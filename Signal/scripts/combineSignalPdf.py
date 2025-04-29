#!/usr/bin/env python
# Combines all signal model PDFs of the specified category without using coefficients

import os, re, sys
import glob
import ROOT
from optparse import OptionParser

def get_options():
    parser = OptionParser()
    parser.add_option("--cat", dest='cat', default='RECO_0J_PTH_0_10_Tag0', help="RECO category for which to combine PDFs")
    parser.add_option("--inputExt", dest='inputExt', default='packaged', help="Folder extension of the input signal model")
    parser.add_option("--outputExt", dest='outputExt', default='combinedPDFs', help="Output folder extension")
    return parser.parse_args()
(opt,args) = get_options()

def rooiter(x):
    iter = x.iterator()
    ret = iter.Next()
    while ret:
        yield ret
        ret = iter.Next()

# Create output workspace
print " --> Creating combined PDF workspace"
combinedWS = ROOT.RooWorkspace("wsig_13TeV", "wsig_13TeV")
combinedWS.imp = getattr(combinedWS, "import")

# Create Higgs mass variable
MH = ROOT.RooRealVar("MH", "MH", 125.0, 120.0, 130.0)
combinedWS.imp(MH)

# List to store all PDFs
all_pdfs = []

# Find all input files
input_files = glob.glob("./outdir_%s/CMS-HGG_sigfit_%s_%s.root" % (opt.inputExt, opt.inputExt, opt.cat))
if not input_files:
    print "Warning: No input files found for category %s" % (opt.cat)
    sys.exit(1)

# Process each input file
for input_file in input_files:
    fin = ROOT.TFile(input_file)
    wsin = fin.Get("wsig_13TeV")
    if not wsin: continue

    mass_point = wsin.var("MH")
    mass_point.setVal(125.0)
    
    # Extract all objects
    allVars = {}
    allPdfs = {}
    
    for _var in rooiter(wsin.allVars()): allVars[_var.GetName()] = _var
    for _pdf in rooiter(wsin.allPdfs()): allPdfs[_pdf.GetName()] = _pdf
    
    # Import all variables
    for varName, var in allVars.iteritems():
        combinedWS.imp(var, ROOT.RooFit.RecycleConflictNodes(), ROOT.RooFit.Silence())
    
    # Find all signal PDFs for this category
    for pdfName, pdf in allPdfs.iteritems():
        if "extenddcb_" in pdfName and opt.cat in pdfName and "ThisLumi" not in pdfName:
            # Store PDF
            all_pdfs.append(pdf)
            print "Found PDF: %s" % pdfName

# Before combining, evaluate each PDF at different mass points
print " --> Evaluating individual PDFs at different mass points (80-180 GeV)"

# Create mass variable if not already in workspace
mgg = combinedWS.var("CMS_hgg_mass")
if not mgg:
    print "Creating CMS_hgg_mass variable"
    mgg = ROOT.RooRealVar("CMS_hgg_mass", "Diphoton mass", 100, 180, "GeV")
    combinedWS.imp(mgg)

# Set mass range for evaluation
min_mass = 80.0
max_mass = 180.0
step_size = 1.0

# Ensure output directory exists
if not os.path.isdir("outdir_%s" % opt.outputExt):
    os.system("mkdir -p outdir_%s" % opt.outputExt)

# Prepare output file for individual PDFs
indiv_values_file = "./outdir_%s/individual_pdf_values_%s.txt" % (opt.outputExt, opt.cat)
with open(indiv_values_file, "w") as fout:
    fout.write("# Mass point (GeV)")
    for i, pdf in enumerate(all_pdfs):
        fout.write("\tPDF_%d" % i)
    fout.write("\n")
    
    # Ensure mgg range covers our needed range
    orig_min = mgg.getMin()
    orig_max = mgg.getMax()
    mgg.setRange(min_mass, max_mass)
    
    # Sample key mass points to show in terminal (to avoid too much output)
    sample_masses = [80.0, 100.0, 120.0, 125.0, 130.0, 150.0, 180.0]
    
    # Print header for terminal output
    print "Mass (GeV)",
    for i in range(len(all_pdfs)):
        print "\tPDF_%d" % i,
    print ""
    print "-" * (15 + 15 * len(all_pdfs))
    
    # Calculate PDF values at each mass point
    for mass in [min_mass + i * step_size for i in range(int((max_mass - min_mass) / step_size) + 1)]:
        mgg.setVal(mass)
        
        # Write to file for all mass points
        fout.write("%.1f" % mass)
        
        # For each PDF, calculate normalized value
        for pdf in all_pdfs:
            # Normalize the PDF
            norm = pdf.createIntegral(ROOT.RooArgSet(mgg))
            norm_val = norm.getVal()
            print "norm_val: %s, valv: %s" % (pdf.getNorm(ROOT.RooArgSet(mgg)), pdf.getValV(ROOT.RooArgSet(mgg)))
            
            # Get PDF value at this mass
            pdf_val = pdf.getVal()
            normalized_val = pdf_val / norm_val if norm_val > 0 else 0
            
            # Write to file
            fout.write("\t%.8e" % normalized_val)
            
            break
        
        # End the line in the file
        fout.write("\n")
    
    # End the last line in terminal output
    print ""
    
    # Restore original range
    mgg.setRange(orig_min, orig_max)

print "Individual PDF values saved to file: %s" % indiv_values_file
print ""

# Create a single combined PDF with all signal components
print " --> Creating combined PDF with %d components" % len(all_pdfs)

# Create RooArgList for PDFs
pdf_list = ROOT.RooArgList()

# Add all PDFs
for pdf in all_pdfs:
    pdf_list.add(pdf)

# Create the combined PDF if we have PDFs
if len(all_pdfs) > 0:
    combined_name = "combinedSigPdf_%s" % (opt.cat)
    
    if len(all_pdfs) == 1:
        # If only one PDF, use it directly
        combined_pdf = all_pdfs[0]
        combinedWS.imp(combined_pdf, ROOT.RooFit.Rename(combined_name))
    else:
        # Create a single combined PDF with all components without coefficients
        # This will create a PDF that is the sum of all components with equal weights
        combined_pdf = ROOT.RooAddPdf(combined_name, "Combined signal PDF (equal weights)", pdf_list)
        combinedWS.imp(combined_pdf)
    
    print "Created combined PDF: %s with %d components" % (combined_name, len(all_pdfs))
else:
    print "Error: No PDFs found for combining"
    sys.exit(1)

# Ensure output directory exists
if not os.path.isdir("outdir_%s" % opt.outputExt):
    os.system("mkdir -p outdir_%s" % opt.outputExt)

# Create a simple validation plot
print " --> Creating validation plot"

# Get observable from workspace
mgg = combinedWS.var("CMS_hgg_mass")
if not mgg:
    print "Warning: Could not find CMS_hgg_mass in workspace, creating it"
    mgg = ROOT.RooRealVar("CMS_hgg_mass", "Diphoton mass", 100, 180, "GeV")
    combinedWS.imp(mgg)

# Create a frame for the plot
frame = mgg.frame(ROOT.RooFit.Title("Combined Signal Model"))

# Plot the combined PDF
combined_pdf.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kBlue))

# Create canvas and draw
c1 = ROOT.TCanvas("c1", "Combined Signal Model", 800, 600)
frame.Draw()

# Add text with info
text = ROOT.TPaveText(0.65, 0.75, 0.89, 0.89, "NDC")
text.SetBorderSize(0)
text.SetFillColor(0)
text.SetTextAlign(12)
text.SetTextFont(42)
text.SetTextSize(0.03)
text.AddText("Category: %s" % opt.cat)
text.AddText("Total components: %d" % len(all_pdfs))
text.AddText("Equal weights")
text.Draw("same")

# Save plot
plot_filename = "./outdir_%s/combined_signal_plot_%s.png" % (opt.outputExt, opt.cat)
c1.SaveAs(plot_filename)
print "Validation plot saved to: %s" % plot_filename

# Save workspace to output file
output_filename = "./outdir_%s/CMS-HGG_combinedPDFs_%s.root" % (opt.outputExt, opt.cat)
f = ROOT.TFile(output_filename, "RECREATE")
combinedWS.Write()
f.Close()

# Print access information
print "Done! Combined PDF saved to output file, accessible via:"
print "  TFile *f = TFile::Open(\"%s\");" % output_filename
print "  RooWorkspace *w = (RooWorkspace*)f->Get(\"wsig_13TeV\");"
print "  RooAbsPdf *pdf = w->pdf(\"%s\");" % combined_name

# Now reopen the saved PDF file and read values from 80-180 GeV in 1 GeV steps
print "\n --> Reading saved PDF file, outputting values in 80-180 GeV range"

# Already saved above, now use directly
min_mass = 80.0
max_mass = 180.0
step_size = 1.0

# Prepare output file
output_values_file = "./outdir_%s/pdf_values_%s.txt" % (opt.outputExt, opt.cat)
with open(output_values_file, "w") as fout:
    fout.write("# Mass point (GeV) \t PDF value\n")
    
    # Ensure mgg range covers our needed range
    orig_min = mgg.getMin()
    orig_max = mgg.getMax()
    mgg.setRange(min_mass, max_mass)
    
    # Normalize PDF
    norm = combined_pdf.createIntegral(ROOT.RooArgSet(mgg))
    norm_val = norm.getVal()
    
    # Calculate PDF values within specified range with given step size
    print "Mass point (GeV) \t PDF value"
    print "-" * 30
    
    for mass in [min_mass + i * step_size for i in range(int((max_mass - min_mass) / step_size) + 1)]:
        mgg.setVal(mass)
        pdf_val = combined_pdf.getVal()
        normalized_val = pdf_val / norm_val if norm_val > 0 else 0
        
        # Print to terminal
        print "%.1f \t\t %.8e" % (mass, normalized_val)
        
        # Write to file
        fout.write("%.1f \t %.8e\n" % (mass, normalized_val))
    
    # Restore original range
    mgg.setRange(orig_min, orig_max)

print "\nPDF values saved to file: %s" % output_values_file
