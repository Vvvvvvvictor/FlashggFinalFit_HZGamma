# Script to calculate photon systematics
# * Run script once per category, loops over signal processes
# * Output is pandas dataframe 

print " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ HGG PHOTON SYST CALCULATOR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "
import ROOT
import pandas as pd
import pickle
import os, sys
from optparse import OptionParser
import glob
import re

# From tools
from plottingTools import * #getEffSigma function
from commonTools import *
from commonObjects import *

def leave():
  print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ HGG PHOTON SYST CALCULATOR (END) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "
  sys.exit(1)

def get_options():
  parser = OptionParser()
  parser.add_option("--xvar", dest='xvar', default='CMS_hgg_mass', help="Observable")
  parser.add_option("--cat", dest='cat', default='', help="RECO category")
  # parser.add_option("--procs", dest='procs', default='', help="Signal processes")
  parser.add_option("--flavs", dest='flavs', default='', help="Flavours")
  parser.add_option("--ext", dest='ext', default='', help="Extension")
  # parser.add_option("--inputWSDir", dest='inputWSDir', default='', help="Input flashgg WS directory")
  parser.add_option("--scales", dest='scales', default='CMS_scale_e,CMS_scale_g,CMS_scale_m', help="Shape systematics: scales")
  # parser.add_option("--scalesCorr", dest='scalesCorr', default='', help='Photon shape systematics: scalesCorr')
  # parser.add_option("--scalesGlobal", dest='scalesGlobal', default='', help='Photon shape systematics: scalesGlobal')
  parser.add_option("--smears", dest='smears', default='CMS_res_e,CMS_res_g,CMS_res_m', help='Shape systematics: smears')
  parser.add_option("--nBins", dest='nBins', default=340, type='int', help='Number of bins in histograms')
  parser.add_option("--thresholdMean", dest='thresholdMean', default=0.05, type='float', help='Reject mean variations if larger than thresholdMean')
  parser.add_option("--thresholdSigma", dest='thresholdSigma', default=0.5, type='float', help='Reject mean variations if larger than thresholdSigma')
  parser.add_option("--thresholdRate", dest='thresholdRate', default=0.05, type='float', help='Reject mean variations if larger than thresholdRate')
  return parser.parse_args()
(opt,args) = get_options()

catMap_sync = {
    "ggH0": "ggf4", "ggH1": "ggf3", "ggH2": "ggf2", "ggH3": "ggf1",
    "VBF0": "vbf4", "VBF1": "vbf3", "VBF2": "vbf2", "VBF3": "vbf1"
}
flavMap_sync = {}

# RooRealVar to fill histograms
xvarName = opt.xvar + catMap_sync[opt.cat]
mgg = ROOT.RooRealVar(xvarName,xvarName,125)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to extact histograms from WS
def getHistograms( _ws, _nominalDataName, _sname ):
  _hists = {}
  # Define histograms
  for htype in ['nominal','Up','Down']:
    if htype == 'nominal': _hists[htype] = ROOT.TH1F(htype,htype,opt.nBins,100,180)
    else: _hists[htype] = ROOT.TH1F("%s%s"%(_sname,htype),"%s_%s"%(_sname,htype),opt.nBins,100,180)
  # Extract nominal RooDataSet and syst RooDataHists
  rds_nominal = _ws.data("%s_nominal"%(_nominalDataName))
  rdh_up = _ws.data("%s_%sUp"%(_nominalDataName,_sname))
  rdh_down = _ws.data("%s_%sDown"%(_nominalDataName,_sname))
  # Check if not NONE type and fill histograms
  if rds_nominal: rds_nominal.fillHistogram(_hists['nominal'],ROOT.RooArgList(mgg))
  else:
    print " --> [ERROR] Could not extract nominal RooDataSet: %s. Leaving"%_nominalDataName
    sys,exit(1)
  if rdh_up: rdh_up.fillHistogram(_hists['Up'],ROOT.RooArgList(mgg))
  else:
    print " --> [ERROR] Could not extract RooDataHist (%s,up) for %s. Leaving"%(_sname,_nominalDataName)
    sys,exit(1)
  if rdh_down: rdh_down.fillHistogram(_hists['Down'],ROOT.RooArgList(mgg))
  else:
    print " --> [ERROR] Could not extract RooDataHist (%s,down) for %s. Leaving"%(_sname,_nominalDataName)
    sys,exit(1) 
  return _hists

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Functions to extract mean, sigma and rate variations
def getMeanVar(_hists):
  mu, muVar = {}, {}
  for htype,h in _hists.iteritems(): mu[htype] = h.GetMean()
  if mu['nominal']==0: return 0
  for htype in ['Up','Down']: muVar[htype] = (mu[htype]-mu['nominal'])/mu['nominal']
  x = (abs(muVar['Up'])+abs(muVar['Down']))/2
  # Check for NaN
  if x!=x: return 0
  else: return min(x,opt.thresholdMean)

def getSigmaVar(_hists):
  sigma, sigmaVar = {}, {}
  for htype,h in _hists.iteritems(): sigma[htype] = getEffSigma(h)
  if sigma['nominal']==0: return 0
  for htype in ['Up','Down']: sigmaVar[htype] = (sigma[htype]-sigma['nominal'])/sigma['nominal']
  x = (abs(sigmaVar['Up'])+abs(sigmaVar['Down']))/2
  if x!=x: return 0
  else: return min(x,opt.thresholdSigma)

def getRateVar(_hists):
  rate, rateVar = {}, {}
  for htype,h in _hists.iteritems(): rate[htype] = h.Integral()
  # Shape variations can both be one sided therefore use midpoint as nominal
  rate['midpoint'] = 0.5*(rate['Up']+rate['Down'])
  if rate['midpoint']==0: return 0
  for htype in ['Up','Down']: rateVar[htype] = (rate[htype]-rate['midpoint'])/rate['midpoint']
  x = (abs(rateVar['Up'])+abs(rateVar['Down']))/2
  if x!=x: return 0
  else: return min(x,opt.thresholdRate)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define dataFrame
columns_data = ['cat','flav','inputWSFile','nominalDataName']
for stype in ['scales','smears']:
  systs = getattr( opt, stype )
  for s in systs.split(","):
    if s == '': continue
    for x in ['mean','sigma','rate']: columns_data.append("%s_%s"%(s,x))
data = pd.DataFrame( columns=columns_data ) 

# Loop over processes and add row to dataframe
# Glob M125 filename
_WSFileName = "/afs/cern.ch/user/m/mioshiro/public/test_datacard11_rawdata.root" # /eos/project/h/htozg-dy-privatemc/mioshiro/test_datacard11_rawdata.root
print("Input ws file: %s"%_WSFileName)
for _flav in opt.flavs.split(","):
  _nominalDataName = "mcdata_Htomm_cat_%s"%(catMap_sync[opt.cat])
  data = data.append({'proc':"all",'cat':opt.cat,'flav':"all",'inputWSFile':_WSFileName,'nominalDataName':_nominalDataName}, ignore_index=True, sort=False)

# Loop over rows in dataFrame and open ws
for ir,r in data.iterrows():

  print " --> Processing (%s,%s)"%(opt.cat,r['flav'])

  # Open ROOT file and extract workspace
  f = ROOT.TFile(r['inputWSFile'])
  _WSName = "WS_Htomm_cat_%s"%(catMap_sync[opt.cat])
  print("Input ws: %s"%_WSName)
  inputWS = f.Get(_WSName)
 
  # Loop over scale and smear systematics
  for stype in ['scales','smears']:
    for s in getattr(opt,stype).split(","):
      if s == '': continue
      # sname = "%s%s"%(inputNuisanceExtMap[stype],s)
      # Following has been modified by JLS just for the early analysis, we need to revise the interfacing
      if s=='CMS_scale_e' and stype=='scales':
        sname = 'CMS_scale_e'
      elif s=='CMS_scale_m' and stype=='scales':
        sname = 'CMS_scale_m'
      elif s=='CMS_scale_g' and stype=='scales':
        sname = 'CMS_scale_g'
      elif s=='CMS_res_e' and stype=='smears':
        sname = 'CMS_res_e'
      elif s=='CMS_res_m' and stype=='smears':
        sname = 'CMS_res_m'
      elif s=='CMS_res_g' and stype=='smears':
        sname = 'CMS_res_g'
      else:
        sname = "%s%s"%(inputNuisanceExtMap[stype],s)
      print "    * Systematic : %s (%s)"%(sname,stype)
      hists = getHistograms(inputWS,r['nominalDataName'],sname) # sname = CMS_{scale/res}_{e/m/g}
      # If nominal yield = 0:
      if hists['nominal'].Integral() == 0: _meanVar, _sigmaVar, _rateVar = 0, 0, 0
      else:
        _meanVar = getMeanVar(hists)
        _sigmaVar = getSigmaVar(hists)
        _rateVar = getRateVar(hists)
      # Add values to dataFrame
      if "scale" in sname: _sigmaVar = 0
      if "res" in sname: _meanVar = 0
      if "_m" in sname and r['flav'] == "ele": _meanVar=_sigmaVar=_rateVar=0
      if "_e" in sname and r['flav'] == "mu": _meanVar=_sigmaVar=_rateVar=0
      data.at[ir,'%s_mean'%(sname)] = _meanVar
      data.at[ir,'%s_sigma'%(sname)] = _sigmaVar
      data.at[ir,'%s_rate'%(sname)] = _rateVar

      # Delete histograms
      for h in hists.itervalues(): h.Delete()

  # Delete ws and close file
  inputWS.Delete()
  f.Close()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Output dataFrame as pickle file to be read in by signalFit.py
if not os.path.isdir("%s/output_%s"%(swd__,opt.ext)): os.system("mkdir %s/output_%s"%(swd__,opt.ext))
if not os.path.isdir("%s/output_%s/calcPhotonSyst_hmm"%(swd__,opt.ext)): os.system("mkdir %s/output_%s/calcPhotonSyst_hmm"%(swd__,opt.ext))
if not os.path.isdir("%s/output_%s/calcPhotonSyst_hmm/pkl"%(swd__,opt.ext)): os.system("mkdir %s/output_%s/calcPhotonSyst_hmm/pkl"%(swd__,opt.ext))
with open("%s/output_%s/calcPhotonSyst_hmm/pkl/%s.pkl"%(swd__,opt.ext,opt.cat),"wb") as f: pickle.dump(data,f) 
with open("%s/output_%s/calcPhotonSyst_hmm/pkl/%s.txt"%(swd__,opt.ext,opt.cat),"wb") as f: f.write(data.to_string(index=False))
print " --> Successfully saved photon systematics as pkl file: %s/output_%s/calcPhotonSyst_hmm/pkl/%s.pkl"%(swd__,opt.ext,opt.cat)
print(data)