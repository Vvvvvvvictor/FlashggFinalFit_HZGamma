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
        if "extendhggpdfsmrel_" in pdfName and opt.cat in pdfName:
            # Store PDF
            all_pdfs.append(pdf)
            print "Found PDF: %s" % pdfName

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
