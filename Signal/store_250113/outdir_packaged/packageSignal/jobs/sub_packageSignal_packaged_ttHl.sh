#!/bin/bash
ulimit -s unlimited
set -e
cd /eos/home-j/jiehan/finalfit_102X/CMSSW_10_2_13/src
export SCRAM_ARCH=slc7_amd64_gcc700
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd /eos/home-j/jiehan/finalfit_102X/CMSSW_10_2_13/src/flashggFinalFit/Signal
export PYTHONPATH=$PYTHONPATH:/eos/home-j/jiehan/finalfit_102X/CMSSW_10_2_13/src/flashggFinalFit/tools:/eos/home-j/jiehan/finalfit_102X/CMSSW_10_2_13/src/flashggFinalFit/Signal/tools

python /eos/home-j/jiehan/finalfit_102X/CMSSW_10_2_13/src/flashggFinalFit/Signal/scripts/packageSignal.py --cat ttHl --outputExt packaged --massPoints 120,125,130 --exts fiducial_2016preVFP,fiducial_2016postVFP,fiducial_2017,fiducial_2018 --mergeYears
