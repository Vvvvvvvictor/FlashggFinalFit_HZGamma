# Script for making signal model plot
import os, sys
import ROOT
import re, glob
import json
from optparse import OptionParser

from commonTools import *
from commonObjects import *
from tools.plottingTools import *

def get_options():
  parser = OptionParser()
  parser.add_option('--procs', dest='procs', default='all', help="Comma separated list of processes to include. all = sum all signal procs")  
  parser.add_option('--years', dest='years', default='2016,2017,2018', help="Comma separated list of years to include")  
  parser.add_option('--cats', dest='cats', default='', help="Comma separated list of analysis categories to include. all = sum of all categories, wall = weighted sum of categories (requires S/S+B from ./Plots/getCatInfo.py)")
  parser.add_option('--flavs', dest='flavs', default='all', help="Comma separated list of flavors to include. all = sum of all flavors")
  parser.add_option('--loadCatWeights', dest='loadCatWeights', default='', help="Load S/S+B weights for analysis categories (path to weights json file)")
  parser.add_option('--ext', dest='ext', default='test', help="Extension: defines output dir where signal models are saved")
  parser.add_option('--outputPath', dest='outputPath', default='.', help="Output path")
  parser.add_option("--xvar", dest="xvar", default='CMS_hgg_mass:m_{ll#gamma}:GeV', help="x-var (name:title:units)")
  parser.add_option("--mass", dest="mass", default='125', help="Mass of datasets")
  parser.add_option("--MH", dest="MH", default='125', help="Higgs mass (for pdf)")
  parser.add_option("--nBins", dest="nBins", default=170, type='int', help="Number of bins")
  parser.add_option("--pdf_nBins", dest="pdf_nBins", default=3200, type='int', help="Number of bins")
  parser.add_option("--threshold", dest="threshold", default=0.001, type='float', help="Threshold to prune process from plot default = 0.1% of total category norm")
  parser.add_option("--translateCats", dest="translateCats", default=None, help="JSON to store cat translations")
  parser.add_option("--translateProcs", dest="translateProcs", default=None, help="JSON to store proc translations")
  parser.add_option("--label", dest="label", default='Simulation Preliminary', help="CMS Sub-label")
  parser.add_option("--doFWHM", dest="doFWHM", default=False, action='store_true', help="Do FWHM")
  return parser.parse_args()
(opt,args) = get_options()

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

# Extract input files: for first file extract xvar
inputFiles = od()
citr = 0
if opt.cats in ['all','wall']:
  fs = glob.glob("%s/%s/outdir_%s/CMS-HGG_sigfit_%s_*.root"%(swd__,opt.outputPath,opt.ext,opt.ext))
  # print(fs)
  for f in fs:
    cat = re.sub(".root","",f.split("/")[-1].split("_%s_"%opt.ext)[-1])
    inputFiles[cat] = f
    if citr == 0:
      w = ROOT.TFile(f).Get("wsig_13TeV")
      xvar = w.var(opt.xvar.split(":")[0])
      xvar.setPlotLabel(opt.xvar.split(":")[1])
      xvar.setUnit(opt.xvar.split(":")[2])
      alist = ROOT.RooArgList(xvar)
    citr += 1
else:
  for cat in opt.cats.split(","):
    f = "%s/%s/outdir_%s/CMS-HGG_sigfit_%s_%s.root"%(swd__,opt.outputPath,opt.ext,opt.ext,cat)
    # print(f)
    inputFiles[cat] = f
    if citr == 0:
      w = ROOT.TFile(f).Get("wsig_13TeV")
      xvar = w.var(opt.xvar.split(":")[0])
      xvar.setPlotLabel(opt.xvar.split(":")[1])
      xvar.setUnit(opt.xvar.split(":")[2])
      alist = ROOT.RooArgList(xvar)
    citr += 1

# Load cat S/S+B weights
if opt.loadCatWeights != '':
  with open( opt.loadCatWeights ) as jsonfile: catsWeights = json.load(jsonfile)

# Define dict to store data histogram and inclusive + per-year pdf histograms
hists = od()
hists['data'] = xvar.createHistogram("h_data", ROOT.RooFit.Binning(opt.nBins))

# Loop over files
for cat,f in inputFiles.iteritems():
  print " --> Processing %s: file = %s"%(cat,f)

  # Define cat weight
  wcat = catsWeights[cat] if opt.loadCatWeights != '' else 1.

  # Open signal workspace
  fin = ROOT.TFile(f)
  w = fin.Get("wsig_13TeV")
  w.var("MH").setVal(float(opt.MH))

  # Extract normalisations
  norms = od()
  data_rwgt = od()
  hpdfs = od()
  for year in opt.years.split(","):
    # Determine which processes and flavors to use
    procs_list = []
    flavs_list = []
    
    allNorms = w.allFunctions().selectByName("*%s*normThisLumi"%year)
    if opt.procs == 'all':
      for norm in rooiter(allNorms):
        proc = norm.GetName().split("_%s"%year)[0].split("dcb_")[1]
        if proc not in procs_list:
          procs_list.append(proc)
    else:
      procs_list = opt.procs.split(",")
    
    if opt.flavs == 'all':
      allNorms = w.allFunctions().selectByName("*%s*%s*normThisLumi"%(year, cat))
      for norm in rooiter(allNorms):
        flav = norm.GetName().split("%s_"%cat)[-1].split("_%s"%sqrts__)[0]
        if flav not in flavs_list:
          flavs_list.append(flav)
    else:
      flavs_list = opt.flavs.split(",")
    
    # Common loop for all cases
    for proc in procs_list:
      for flav in flavs_list:
        k = "dcb_%s_%s_%s"%(proc, year, flav)
        _id = "dcb_%s_%s_%s_%s_%s"%(proc, year, cat, flav, sqrts__)
        print("k = %s, _id = %s"%(k, _id))
        norm_func = w.function("%s_normThisLumi"%(_id))
        if norm_func:  # Check if function exists
          norms[k] = norm_func

  # Iterate over norms: extract total category norm
  catNorm = 0
  for k, norm in norms.iteritems():
    proc = k.split("_20")[0].split("dcb_")[1]
    year, flav = k.split("%s_"%proc)[-1].split("_")
    _id = "dcb_%s_%s_%s_%s_%s"%(proc, year, cat, flav, sqrts__)
    w.var("IntLumi").setVal(lumiScaleFactor*lumiMap[year])
    catNorm += norm.getVal()
    print("k = %s, norm = %f" % (k,norm.getVal()))
  print("catNorm = %f" % catNorm)
  if catNorm == 0: continue

  # Iterate over norms and extract data sets + pdfs
  for k, norm in norms.iteritems():
    proc = k.split("_20")[0].split("dcb_")[1]
    year, flav = k.split("%s_"%proc)[-1].split("_")
    _id = "dcb_%s_%s_%s_%s_%s"%(proc, year, cat, flav, sqrts__)
    w.var("IntLumi").setVal(lumiScaleFactor*lumiMap[year])
    # print("k = %s, proc = %s, year = %s, flav = %s, _id = %s"%(k, proc, year, flav, _id))

    # Prune
    nval = norm.getVal()
    # print("nval = %f" % nval)
    if nval < opt.threshold*catNorm: continue # Prune processes which contribute less that threshold of signal mod

    # Make empty copy of dataset
    d = w.data("sig_mass_m%s_%s"%(opt.mass,_id.split("dcb_")[-1]))
    d_rwgt = d.emptyClone(_id)
    
    # Calc norm factor
    if d.sumEntries() == 0: nf = 0
    else: nf = nval/d.sumEntries()
    # print("nf = %f, nval = %f, d.sumEntries() = %f" % (nf,nval,d.sumEntries()))

    # Fill dataset with correct normalisation + reweight if using cat weights
    for i in range(d.numEntries()):
      p = d.get(i)
      rw, rwe = d.weight()*nf*wcat, d.weightError()*nf*wcat
      d_rwgt.add(p,rw,rwe)
    # Add dataset to container
    data_rwgt[_id] = d_rwgt

    # Extract pdf and create histogram
    pdf = w.pdf("extend%sThisLumi"%(_id)) 
    # print("extend%sThisLumi"%(_id))
    hpdfs[_id] = pdf.createHistogram("h_pdf_%s"%_id,xvar,ROOT.RooFit.Binning(opt.pdf_nBins))
    # print("wcat = %f, float(opt.nBins) = %f" % (wcat,float(opt.nBins)))
    # hpdfs[_id].Scale(wcat*float(opt.nBins)/float(opt.pdf_nBins)) # FIXME: hardcoded 320
    hpdfs[_id].Scale(wcat*float(opt.nBins)/340) # FIXME: hardcoded 340

  # Fill total histograms: data, per-year pdfs and pdfs
  for _id,d in data_rwgt.iteritems(): d.fillHistogram(hists['data'],alist)
  print("data yield: %f" % hists['data'].Integral())

  # Sum pdf histograms
  for _id,p in hpdfs.iteritems():
    print("hpdfs.iteritems(): _id = %s" % (_id))
    if 'pdf' not in hists: 
      hists['pdf'] = p.Clone("h_pdf")
      hists['pdf'].Reset()
    # Fill
    hists['pdf'] += p
  print("yield total: %f " % (hists['pdf'].Integral()*float(opt.nBins)/float(opt.pdf_nBins)))

  # Per-year pdf histograms
  if len(opt.years.split(",")) > 1:
    for year in opt.years.split(","):
      # print(year)
      if 'pdf_%s'%year not in hists:
        hists['pdf_%s'%year] = hists['pdf'].Clone()
        hists['pdf_%s'%year].Reset()
      # Fill
      for _id,p in hpdfs.iteritems():
        if year in _id: 
          # print("_id = %s" % _id)
          hists['pdf_%s'%year] += p
      # print("yield %s: %f" % (year,hists['pdf_%s' % year].Integral()))
   
  # Garbage removal
  for d in data_rwgt.itervalues(): d.Delete()
  for p in hpdfs.itervalues(): p.Delete()
  w.Delete()
  fin.Close()

# Make plot
if not os.path.isdir("%s/%s/outdir_%s/Plots"%(swd__,opt.outputPath,opt.ext)): os.system("mkdir %s/%s/outdir_%s/Plots"%(swd__,opt.outputPath,opt.ext))
plotSignalModel(hists,opt,_outdir="%s/%s/outdir_%s/Plots"%(swd__,opt.outputPath,opt.ext))