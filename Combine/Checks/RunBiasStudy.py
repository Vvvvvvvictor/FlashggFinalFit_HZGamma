#!/usr/bin/env python
# -*- coding: utf-8 -*-

## sample:
# python /afs/cern.ch/user/s/shsong/CMSSW_10_6_20/src/flashggFinalFit/RunbiasStudy.py -d /afs/cern.ch/user/s/shsong/CMSSW_10_6_20/src/flashggFinalFit/Datacard/combine_run2/condor_input/Datacard_M1000_run2allcat.root  -t

# python /afs/cern.ch/user/s/shsong/CMSSW_10_6_20/src/flashggFinalFit/RunbiasStudy.py -d /afs/cern.ch/user/s/shsong/CMSSW_10_6_20/src/flashggFinalFit/Datacard/combine_run2/condor_input/Datacard_M1000_run2allcat.root -f -c "--cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints   --rMin -5 --rMax 5  --freezeParameters MH"
# python /afs/cern.ch/user/s/shsong/CMSSW_10_6_20/src/flashggFinalFit/RunbiasStudy.py -d /afs/cern.ch/user/s/shsong/CMSSW_10_6_20/src/flashggFinalFit/Datacard/combine_run2/condor_input/Datacard_M1000_run2allcat.root -p --gaussianFit

from biasUtils import *
import matplotlib
from array import array
from os import *
import matplotlib.pyplot as plt
matplotlib.use('Agg')
import json

from pdb import set_trace as bt

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-d","--datacard",default="Datacard.root")
parser.add_option("-w","--workspace",default="w")
parser.add_option("-t","--toys",action="store_true", default=False)
parser.add_option("-n","--nToys",default=1000,type="int")
parser.add_option("-f","--fits",action="store_true", default=False)
parser.add_option("-p","--plots",action="store_true", default=False)
parser.add_option("-j","--jobName",default="fiducial",help="Job name for output directories")
parser.add_option("-e","--expectSignal",default='0,1,2,5,10')
parser.add_option("-m","--mH",default=125.,type="float")
parser.add_option("-c","--combineOptions",default="")
parser.add_option("-s","--seed",default=12345,type="int")
parser.add_option("--dryRun",action="store_true", default=False)
parser.add_option("--condor",action="store_true", default=False)
parser.add_option("--poi",default="r")
parser.add_option("--split",default=1000,type="int")
parser.add_option("--selectFunction",default=None)
parser.add_option("--gaussianFit",action="store_true", default=False)
(opts,args) = parser.parse_args()
print
if opts.nToys>opts.split and not opts.nToys%opts.split==0: raise RuntimeError('The number of toys %g needs to be smaller than or divisible by the split number %g'%(opts.nToys, opts.split))

import ROOT as r
r.gROOT.SetBatch(True)
r.gStyle.SetOptStat(2211)

sname=opts.datacard.split('.')[0]

ws = r.TFile(opts.datacard).Get(opts.workspace)

pdfs = rooArgSetToList(ws.allPdfs())
multipdfName = None
for pdf in pdfs:
    if pdf.InheritsFrom("RooMultiPdf"):
        if multipdfName is not None: raiseMultiError() 
        multipdfName = pdf.GetName()
        print 'Conduct bias study for multipdf called %s'%multipdfName
multipdf = ws.pdf(multipdfName)
print

varlist = rooArgSetToList(ws.allCats())
indexName = None
for var in varlist:
    if var.GetName().startswith('pdfindex'):
        if indexName is not None: raiseMultiError()
        indexName = var.GetName()
        print 'Found index called %s'%indexName
print

from collections import OrderedDict as od
indexNameMap = od()
for ipdf in range(multipdf.getNumPdfs()):
    if opts.selectFunction is not None:
        if not multipdf.getPdf(ipdf).GetName().count(opts.selectFunction): continue
    indexNameMap[ipdf] = multipdf.getPdf(ipdf).GetName()

if opts.toys:
    if opts.jobName and not path.isdir('BiasStudy/%s' % opts.jobName):
        system('mkdir -p BiasStudy/%s' % opts.jobName)
    toysDir = 'BiasStudy/%s/BiasToys' % (opts.jobName if opts.jobName else '')
    if not path.isdir(toysDir): system('mkdir -p %s' % toysDir)
    if opts.condor:
        cmds = []
    for sig in opts.expectSignal.split(','):
        toyCmdBase = 'combine -m %.4f -d %s -M GenerateOnly  --toysNoSystematics --expectSignal %.4f -s %g --saveToys %s '%(opts.mH, opts.datacard, float(sig), opts.seed, opts.combineOptions)
        for ipdf,pdfName in indexNameMap.iteritems():
            name = '%s_%g' %(shortName(pdfName), int(sig))
            print("name: ", toyName(name, jobName=opts.jobName))
            if opts.nToys > opts.split:
                for isplit in range(opts.nToys//opts.split):
                    toyCmd = toyCmdBase + ' -t %g -n _%s_split%g --setParameters %s=%g --freezeParameters %s'%(opts.split, name,  isplit, indexName, ipdf, indexName) + ";echo ;" + 'mv higgsCombine_%s* %s'%(name, toyName(name, split=isplit, jobName=opts.jobName))
            else: 
                toyCmd = toyCmdBase + ' -t %g -n _%s --setParameters %s=%g --freezeParameters %s'%(opts.nToys, name, indexName, ipdf, indexName) + ";echo ;" + 'mv higgsCombine_%s* %s'%(name, toyName(name, jobName=opts.jobName))
            if opts.condor:
                cmds.append(toyCmd)
            else:        
                run(toyCmd, dry=opts.dryRun)
            print "toy command line: ", toyCmd
    if opts.condor:
        writeCondorScript(cmds, opts.jobName, dry=opts.dryRun)
print

if opts.fits:
    if opts.jobName and not path.isdir('BiasStudy/%s' % opts.jobName):
        system('mkdir -p BiasStudy/%s' % opts.jobName)
    fitsDir = 'BiasStudy/%s/BiasFits' % (opts.jobName if opts.jobName else '')
    if not path.isdir(fitsDir): system('mkdir -p %s' % fitsDir)
    if opts.condor:
        cmds = []
    for sig in opts.expectSignal.split(','):
        fitCmdBase = 'combine -m %.4f -d %s -M MultiDimFit -P %s --algo singles %s '%(opts.mH, opts.datacard, opts.poi, opts.combineOptions)
        for ipdf,pdfName in indexNameMap.iteritems():
            name = '%s_%g' %(shortName(pdfName), int(sig))
            if opts.nToys > opts.split:
                for isplit in range(opts.nToys//opts.split):
                    fitCmd = fitCmdBase + ' -t %g -n _%s_split%g --toysFile=%s'%(opts.split, name, isplit, toyName(name, split=isplit, jobName=opts.jobName)) + ';echo ;' + 'mv higgsCombine_%s* %s'%(name, fitName(name, split=isplit, jobName=opts.jobName))
                    run(fitCmd, dry=opts.dryRun)
                run('hadd %s BiasFits/*%s*split*.root'%(fitName(name, jobName=opts.jobName), name), dry=opts.dryRun)
            else:
                fitCmd = fitCmdBase + ' -t %g -n _%s --toysFile=%s'%(opts.nToys, name, toyName(name, jobName=opts.jobName)) + ';echo ;' + 'mv higgsCombine_%s* %s'%(name, fitName(name, jobName=opts.jobName))
                run(fitCmd, dry=opts.dryRun)
            if opts.condor:
                cmds.append(fitCmd)
            else:
                run(fitCmd, dry=opts.dryRun)
            print "fit command line: ", fitCmd
    if opts.condor:
        writeCondorScript(cmds, opts.jobName, dry=opts.dryRun, type='fits')

if opts.plots:
    if opts.jobName and not path.isdir('BiasStudy/%s' % opts.jobName):
        system('mkdir -p BiasStudy/%s' % opts.jobName)
    plotsDir = 'BiasStudy/%s/BiasPlots' % (opts.jobName if opts.jobName else '')
    print "!!!!!!!!plotsDir: ", plotsDir
    if not path.isdir(plotsDir): system('mkdir -p %s' % plotsDir)
    pdfnames = []
    means = []
    mean_errors = []
    sigmas = []
    expectSignals = []
    bias_data = {}

    for sig in opts.expectSignal.split(','):
        for ipdf, pdfName in indexNameMap.iteritems():
            name = '%s_%g' %(shortName(pdfName), int(sig))
            tfile = r.TFile(fitName(name, jobName=opts.jobName))
            tree = tfile.Get('limit')
            pullHist = r.TH1F('pullsForTruth_%s'%name, 'Pull distribution using the envelope to fit %s'%name, 80, -4., 4.)
            pullHist.GetXaxis().SetTitle('Pull')
            pullHist.GetYaxis().SetTitle('Entries')
            plotentry = []
            
            for itoy in range(opts.nToys):
                tree.GetEntry(3*itoy)
                if not getattr(tree, 'quantileExpected') == -1: 
                    raiseFailError(itoy, True) 
                    continue
                bf = getattr(tree, 'r')
                tree.GetEntry(3*itoy+1)
                if not abs(getattr(tree, 'quantileExpected') - -0.32) < 0.001: 
                    raiseFailError(itoy, True) 
                    continue
                lo = getattr(tree, 'r')
                tree.GetEntry(3*itoy+2)
                if not abs(getattr(tree, 'quantileExpected') - 0.32) < 0.001: 
                    raiseFailError(itoy, True) 
                    continue
                hi = getattr(tree, 'r')
                diff = bf - int(sig)
                unc = 0.5 * (hi-lo)
                if unc > 0.: 
                    pullHist.Fill(diff/unc)
                    plotentry.append(diff/unc)
                    
            sorted_entry = sorted(plotentry)
            n = len(sorted_entry)
            if n % 2 == 1:
                median = sorted_entry[n // 2]
            else:
                middle1 = sorted_entry[n // 2 - 1]
                middle2 = sorted_entry[n // 2]
                median = (middle1 + middle2) / 2.0
                
            canv = r.TCanvas()
            pullHist.Draw()
            
            mean = 0
            mean_error = 0
            sigma = 0
            if opts.gaussianFit:
                r.gStyle.SetOptFit(111)
                pullHist.Fit('gaus')
                fit_result = pullHist.GetFunction("gaus")
                mean = fit_result.GetParameter(1)
                mean_error = fit_result.GetParError(1)  # Get the fitting error of mean
                sigma = fit_result.GetParameter(2)
            else:
                mean = pullHist.GetMean()
                mean_error = pullHist.GetMeanError()  # Get the statistical error of mean
                sigma = pullHist.GetRMS()
                
            pdfnames.append(pdfName)
            means.append(mean)
            mean_errors.append(mean_error)
            sigmas.append(sigma)
            expectSignals.append(int(sig))
            
            # Save mean, mean_error and sigma data for each function
            bias_data[pdfName] = {
                "mean": mean,
                "mean_error": mean_error,
                "sigma": sigma,
                "expectSignal": int(sig)
            }
            
            canv.SaveAs('%s.pdf' % plotName(name, jobName=opts.jobName))
            canv.SaveAs('%s.png' % plotName(name, jobName=opts.jobName))
    
    # # Save data to JSON file
    # json_output = path.join(plotsDir, 'bias_results.json')  # Use path.join for safer path handling
    # print "Bias results will be saved to: ", json_output
    # # bt()
    # json_file = os.open(json_output, os.O_WRONLY | os.O_CREAT | os.O_TRUNC)
    # json.dump(bias_data, json_file, indent=4)
    # os.close(json_file)
    # print("Bias results saved to: %s" % json_output)
    
    # Create a new plot with mean on X-axis and expectSignal on Y-axis
    canvas = r.TCanvas("canvas", "Bias Study", 800, 600)
    canvas.SetGrid()
    
    # Set appropriate margins
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.05)
    canvas.SetTopMargin(0.1)
    canvas.SetBottomMargin(0.15)
    
    # Create a multi-graph to hold all individual graphs
    mg = r.TMultiGraph()
    mg.SetTitle("Bias Study Results")
    
    # Colors to use for different functions
    colors = [r.kRed, r.kBlue, r.kGreen+2, r.kMagenta, r.kOrange+7, 
             r.kCyan+2, r.kViolet-1, r.kSpring+5, r.kTeal+1, r.kYellow+2]
    
    # Marker styles to use
    markers = [20, 21, 22, 23, 29, 33, 34, 47, 43, 39]
    
    # Create graphs for each function
    graphs = []
    
    for i, (pdfname, mean, mean_error, sigma, expectSignal) in enumerate(zip(pdfnames, means, mean_errors, sigmas, expectSignals)):
        # Create graph for this function
        graph = r.TGraphErrors(1)
        
        # Create a dictionary to map pdfnames to their color and marker indices
        pdfname_to_indices = {}
        unique_pdfnames = list(set(pdfnames))

        for i, unique_pdf in enumerate(unique_pdfnames):
            # Assign a consistent color and marker index to each unique PDF
            pdfname_to_indices[unique_pdf] = {
                'color_idx': i % len(colors),
                'marker_idx': i % len(markers)
            }

        # Get the color and marker index for this PDF
        color_idx = pdfname_to_indices[pdfname]['color_idx']
        marker_idx = pdfname_to_indices[pdfname]['marker_idx']
        
        # Set point
        graph.SetPoint(0, mean, expectSignal)
        graph.SetPointError(0, mean_error, 0)  # x_err=mean_error, y_err=0
        
        # Set graph properties
        graph.SetMarkerStyle(markers[marker_idx])
        graph.SetMarkerSize(1.2)
        graph.SetMarkerColor(colors[color_idx])
        graph.SetLineColor(colors[color_idx])
        graph.SetLineWidth(2)
        
        # Store reference to graph
        graphs.append(graph)
        
        # Add to multigraph
        mg.Add(graph)

    # Calculate appropriate x-axis range
    max_x = max(max([abs(m) for m in means]) * 1.2, 0.2)
    
    # Dynamically adjust y-axis range based on number of legend items
    legend_items = len(unique_pdfnames) + 2 # Number of functions + 2 color region descriptions
    legend_height = legend_items * 0.06  # Each item takes more space    

    # Calculate appropriate y-axis range
    min_y = min(expectSignals) - (max(expectSignals) - min(expectSignals) + 1) * 0.1  # Default lower limit
    max_y = max(expectSignals) + (max(expectSignals) - min(expectSignals) + 1) * (0.2 + 1.1 * legend_height / (1 - legend_height))  # Default upper limit
    
    # Adjust y-axis minimum to ensure enough space for legend
    legend_position_top = 0.9  # Legend top position (NDC coordinates)
    legend_position_bottom = legend_position_top - legend_height
    
    # If legend is too large, adjust y-axis minimum
    if legend_position_bottom < 0.2:  # Ensure bottom margin
        min_y = min(expectSignals) * 0.5  # Further reduce y-axis minimum
    
    # Create an empty histogram as coordinate system frame
    frame = r.TH2F("frame", "Bias Study Results", 100, -max_x, max_x, 100, min_y, max_y)
    frame.SetStats(0)  # Turn off stats box
    
    # Set axis titles and properties
    frame.GetXaxis().SetTitle("(#mu - #tilde{#mu})/#sigma")
    frame.GetYaxis().SetTitle("Expected Signal (#mu)")
    
    # Increase title font size
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetLabelSize(0.045)
    frame.GetYaxis().SetLabelSize(0.045)
    frame.GetXaxis().SetTitleOffset(1.1)
    frame.GetYaxis().SetTitleOffset(1.1)
    
    # Draw frame
    frame.Draw()
    
    # Get current plotting area y-axis range
    y_min = frame.GetYaxis().GetXmin()
    y_max = frame.GetYaxis().GetXmax()
    
    # Add color regions
    # Light yellow region (-0.2, 0.2)
    yellow_box = r.TBox(-0.2, y_min, 0.2, y_max)
    yellow_box.SetFillColorAlpha(r.kYellow, 0.3)
    yellow_box.Draw("same")
    
    # Light green region (-0.14, 0.14)
    green_box = r.TBox(-0.14, y_min, 0.14, y_max)
    green_box.SetFillColorAlpha(r.kGreen, 0.3)
    green_box.Draw("same")
    
    # Set canvas to transparent mode
    canvas.SetFillStyle(0)
    
    # Draw each graph individually instead of using TMultiGraph
    for graph in graphs:
        graph.SetFillStyle(0)  # Set graph fill to transparent
        graph.Draw("P same")  # Draw each graph separately
    
    # Create legend
    # Dynamically adjust legend position and size based on number of items
    legend_width = 0.4
    legend = r.TLegend(0.15, legend_position_bottom, 0.15 + legend_width, legend_position_top)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.04)  # Increase legend text size
    
    # Add entries for each function
    # Create a map of unique PDFs to their corresponding graphs
    pdf_graph_map = {}
    for i, (graph, pdfname) in enumerate(zip(graphs, pdfnames)):
        if pdfname not in pdf_graph_map:
            pdf_graph_map[pdfname] = graph
    
    # Add only unique PDFs to the legend
    for pdfname, graph in pdf_graph_map.items():
        # Use short name for better display
        short_name = shortName(pdfname)
        legend.AddEntry(graph, short_name, "p")
        
    # Add entries for color regions
    legend.AddEntry(yellow_box, "Bias < 0.2", "f")
    legend.AddEntry(green_box, "Bias < 0.14", "f")
    legend.Draw("same")
    
    # Add CMS Preliminary label
    cms_text = r.TLatex()
    cms_text.SetNDC()
    cms_text.SetTextFont(61)
    cms_text.SetTextSize(0.06)  # Larger text
    cms_text.DrawLatex(0.15, 0.92, "CMS")
    
    prelim_text = r.TLatex()
    prelim_text.SetNDC()
    prelim_text.SetTextFont(52)
    prelim_text.SetTextSize(0.045)  # Larger text
    prelim_text.DrawLatex(0.26, 0.92, "Preliminary")
    
    # Add bias value text
    bias_text = r.TLatex()
    bias_text.SetNDC()
    bias_text.SetTextFont(42)
    bias_text.SetTextSize(0.035)  # Larger text
    
    # Calculate text box position and size
    text_x = 0.65
    text_y_top = 0.85
    text_y_step = 0.045  # Larger spacing
    
    for i, (pdfname, mean, mean_error) in enumerate(zip(pdfnames, means, mean_errors)):
        short_name = shortName(pdfname)
        bias_text.DrawLatex(text_x, text_y_top-i*text_y_step, "%s: %.3f #pm %.3f" % (short_name, mean, mean_error))
    
    # Save graph
    canvas.Update()
    output_name = '%s/bias_summary' % plotsDir
    canvas.SaveAs('%s.png' % output_name)
    canvas.SaveAs('%s.pdf' % output_name)
    
    print("Bias summary plot saved to: %s.png/pdf" % output_name)
