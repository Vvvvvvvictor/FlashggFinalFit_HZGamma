#!/bin/bash
ulimit -s unlimited
set -e
cd /eos/home-j/jiehan/finalfit_102X/CMSSW_10_2_13/src
export SCRAM_ARCH=slc7_amd64_gcc700
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd /eos/home-j/jiehan/finalfit_102X/CMSSW_10_2_13/src/flashggFinalFit/Signal
export PYTHONPATH=$PYTHONPATH:/eos/home-j/jiehan/finalfit_102X/CMSSW_10_2_13/src/flashggFinalFit/tools:/eos/home-j/jiehan/finalfit_102X/CMSSW_10_2_13/src/flashggFinalFit/Signal/tools

python /eos/home-j/jiehan/finalfit_102X/CMSSW_10_2_13/src/flashggFinalFit/Signal/scripts/signalFit.py --inputWSDir /eos/home-j/jiehan/root/input_finalfit/signal_2017 --ext fiducial_2017 --proc VBF --cat VBF2 --year 2017 --analysis fiducialAnalysis --massPoints 120,125,130 --scales 'Scale,MuonPt' --scalesCorr 'Material,FNUF' --scalesGlobal '' --smears 'Smearing' --doPlots --beamspotWidthData 3.5 --beamspotWidthMC 3.7

python /eos/home-j/jiehan/finalfit_102X/CMSSW_10_2_13/src/flashggFinalFit/Signal/scripts/signalFit.py --inputWSDir /eos/home-j/jiehan/root/input_finalfit/signal_2017 --ext fiducial_2017 --proc WminusH --cat VBF2 --year 2017 --analysis fiducialAnalysis --massPoints 120,125,130 --scales 'Scale,MuonPt' --scalesCorr 'Material,FNUF' --scalesGlobal '' --smears 'Smearing' --doPlots --beamspotWidthData 3.5 --beamspotWidthMC 3.7

python /eos/home-j/jiehan/finalfit_102X/CMSSW_10_2_13/src/flashggFinalFit/Signal/scripts/signalFit.py --inputWSDir /eos/home-j/jiehan/root/input_finalfit/signal_2017 --ext fiducial_2017 --proc WplusH --cat VBF2 --year 2017 --analysis fiducialAnalysis --massPoints 120,125,130 --scales 'Scale,MuonPt' --scalesCorr 'Material,FNUF' --scalesGlobal '' --smears 'Smearing' --doPlots --beamspotWidthData 3.5 --beamspotWidthMC 3.7

python /eos/home-j/jiehan/finalfit_102X/CMSSW_10_2_13/src/flashggFinalFit/Signal/scripts/signalFit.py --inputWSDir /eos/home-j/jiehan/root/input_finalfit/signal_2017 --ext fiducial_2017 --proc ZH --cat VBF2 --year 2017 --analysis fiducialAnalysis --massPoints 120,125,130 --scales 'Scale,MuonPt' --scalesCorr 'Material,FNUF' --scalesGlobal '' --smears 'Smearing' --doPlots --beamspotWidthData 3.5 --beamspotWidthMC 3.7

python /eos/home-j/jiehan/finalfit_102X/CMSSW_10_2_13/src/flashggFinalFit/Signal/scripts/signalFit.py --inputWSDir /eos/home-j/jiehan/root/input_finalfit/signal_2017 --ext fiducial_2017 --proc ggH --cat VBF2 --year 2017 --analysis fiducialAnalysis --massPoints 120,125,130 --scales 'Scale,MuonPt' --scalesCorr 'Material,FNUF' --scalesGlobal '' --smears 'Smearing' --doPlots --beamspotWidthData 3.5 --beamspotWidthMC 3.7

python /eos/home-j/jiehan/finalfit_102X/CMSSW_10_2_13/src/flashggFinalFit/Signal/scripts/signalFit.py --inputWSDir /eos/home-j/jiehan/root/input_finalfit/signal_2017 --ext fiducial_2017 --proc ttH --cat VBF2 --year 2017 --analysis fiducialAnalysis --massPoints 120,125,130 --scales 'Scale,MuonPt' --scalesCorr 'Material,FNUF' --scalesGlobal '' --smears 'Smearing' --doPlots --beamspotWidthData 3.5 --beamspotWidthMC 3.7

