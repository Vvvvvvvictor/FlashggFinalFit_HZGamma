#!/bin/bash
ulimit -s unlimited
set -e
cd /eos/home-j/jiehan/finalfit_102X/CMSSW_10_2_13/src
export SCRAM_ARCH=slc7_amd64_gcc700
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd /eos/home-j/jiehan/finalfit_102X/CMSSW_10_2_13/src/flashggFinalFit/Datacard
export PYTHONPATH=$PYTHONPATH:/eos/home-j/jiehan/finalfit_102X/CMSSW_10_2_13/src/flashggFinalFit/tools:/eos/home-j/jiehan/finalfit_102X/CMSSW_10_2_13/src/flashggFinalFit/Datacard/tools

python /eos/home-j/jiehan/finalfit_102X/CMSSW_10_2_13/src/flashggFinalFit/Datacard/makeYields.py --cat VHlep --procs auto --ext others --mass 125 --inputWSDirMap 2016preVFP=/eos/home-j/jiehan/root/input_finalfit/signal_2016preVFP/,2016postVFP=/eos/home-j/jiehan/root/input_finalfit/signal_2016postVFP/,2017=/eos/home-j/jiehan/root/input_finalfit/signal_2017/,2018=/eos/home-j/jiehan/root/input_finalfit/signal_2018/ --sigModelWSDir ./Models/signal --sigModelExt packaged --bkgModelWSDir ./Models/background --bkgModelExt multipdf  --mergeYears --skipZeroes --doSystematics
