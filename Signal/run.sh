#!/bin/bash

clear
cd ../
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsenv
echo $(pwd)
source $(pwd)/setup.sh
mainPath=$(pwd)
inputExt="25Dec"
inputpath="input_${inputExt}"
outputExt="25Dec"
outputpath="output_${outputExt}"

BackgroundWSPath="${mainPath}/Inputdata/${inputpath}/background"
BackgroundNtuplePath="${mainPath}/Inputdata/${inputpath}/fitting_bkg"
SignalNtuplePath="${mainPath}/Inputdata/${inputpath}/fitting_signal"
SignalWSPath="${mainPath}/Inputdata/${inputpath}/signal"

# SignalProcs=("ggH" "VBF" "WH" "ZH" "ttH") # "ggH" "VBF" "WplusH" "WminusH" "ZH" "ttH"
SignalProcs=("ggH" "VBF") # "ggH" "VBF"
mass_points=("120" "125" "130")
years=("2016preVFP" "2016postVFP" "2017" "2018" "2022preEE" "2022postEE" "2023preBPix" "2023postBPix") #"2016preVFP" "2016postVFP" "2017" "2018" "2022preEE" "2022postEE" "2023preBPix" "2023postBPix"
# years=("all") #"2016preVFP" "2016postVFP" "2017" "2018" "2022preEE" "2022postEE" "2023preBPix" "2023postBPix" "all"

#########################################
# Signal
#########################################
cd ${mainPath}/Signal
# fTest # no need to do when using dcb model
# for year in "${years[@]}"; do
#     python RunSignalScripts.py --inputConfig config_${year}_fiducial.py --mode fTest --modeOpts "--doPlots"
# done
# # syst
# for year in "${years[@]}"; do
#     python RunSignalScripts.py --inputConfig config_${year}_fiducial.py --mode calcPhotonSyst
# done

# signalfit  ## InitValue(Range) can be set in tools/simultaneousFit.py#L10(L73)
# for year in "${years[@]}"; do
#     echo "Running signal fitting for year $year"
#     # sed -i "s#input_jiehan#input_final#g" config_${year}_fiducial.py
#     python RunSignalScripts.py --inputConfig config_${year}_fiducial.py --mode signalFit --groupSignalFitJobsByCat --outputPath ${outputpath} --modeOpts "--doPlots --beamspotWidthData 3.5 --beamspotWidthMC 3.7 --useDCB --skipSystematics"
# done

# # packaged
# echo "Packaging signal models"
# python RunPackager.py --cats ggH0,ggH1,ggH2,ggH3,VBF0,VBF1,VBF2,VBF3 --exts fiducial_2016preVFP,fiducial_2016postVFP,fiducial_2017,fiducial_2018,fiducial_2022preEE,fiducial_2022postEE,fiducial_2023preBPix,fiducial_2023postBPix --mergeYears --batch local --outputExt packaged --outputPath ${outputpath}
# ggH0,ggH1,ggH2,ggH3,VBF0,VBF1,VBF2,VBF3,VHlep,ZHinv,ttHl,ttHh
# Incl0,Incl1,Incl2,Incl3,Incl4,Incl5,Incl6,Incl7

# # signal model plotting
# echo "Plotting signal models"
# python RunPlotter.py --procs all --cats all --years 2016preVFP,2016postVFP,2017,2018,2022preEE,2022postEE,2023preBPix,2023postBPix --ext packaged --outputPath ${outputpath}

# for proc in ggH VBF; do
#     echo "Plotting signal models for $proc production mode"
#     python RunPlotter.py --procs $proc --cats all --years 2016preVFP,2016postVFP,2017,2018,2022preEE,2022postEE,2023preBPix,2023postBPix --ext packaged --outputPath ${outputpath}
# done

# for cat in ggH0 ggH1 ggH2 ggH3 VBF0 VBF1 VBF2 VBF3; do
#     echo "Plotting signal models for $cat category"
#     python RunPlotter.py --procs all --cats $cat --years 2016preVFP,2016postVFP,2017,2018,2022preEE,2022postEE,2023preBPix,2023postBPix --ext packaged --outputPath ${outputpath}
# done

# for flav in ele mu; do
#     echo "Plotting signal models for $flav lepton flavor"
#     python RunPlotter.py --procs all --cats all --flavs $flav --years 2016preVFP,2016postVFP,2017,2018,2022preEE,2022postEE,2023preBPix,2023postBPix --ext packaged --outputPath ${outputpath}
# done

### Merge procs/years/flavs     
### Merge all   ## For hmm fitting
# for cat in ggH0 ggH1 ggH2 ggH3 VBF0 VBF1 VBF2 VBF3; do
#     python scripts/calcPhotonSyst_hmm_sync.py --xvar mllg_cat_ --cat ${cat} --ext ${outputExt}_hmm --scales CMS_scale_e,CMS_scale_g,CMS_scale_m --smears CMS_res_e,CMS_res_g,CMS_res_m
# done
# python RunSignalScripts.py --inputConfig config_all_fiducial.py --mode signalFit --groupSignalFitJobsByCat --outputPath ${outputpath}_hmm --modeOpts "--ext ${outputExt}_hmm --doPlots --skipBeamspotReweigh --useDCB" #  --skipSystematics
# python RunPackager.py --cats ggH0,ggH1,ggH2,ggH3,VBF0,VBF1,VBF2,VBF3 --exts fiducial_all --mergeYears --batch local --outputExt packaged --outputPath ${outputpath}_hmm
# python RunPlotter.py --procs all --cats all --years all --ext packaged --outputPath ${outputpath}_hmm
# for cat in ggH0 ggH1 ggH2 ggH3 VBF0 VBF1 VBF2 VBF3; do
#     echo "Plotting signal models for $cat category"
#     python RunPlotter.py --procs all --cats $cat --years all --ext packaged_all --outputPath ${outputpath}_hmm
# done

### Merge proc
# for year in "${years[@]}"; do
#     echo "Running signal fitting for year $year"
#     sed -i "s#'procs':'auto',#'procs':'all',#g" config_${year}_fiducial.py
#     python RunSignalScripts.py --inputConfig config_${year}_fiducial.py --mode signalFit --groupSignalFitJobsByCat --outputPath ${outputpath}_proc --modeOpts "--doPlots --beamspotWidthData 3.5 --beamspotWidthMC 3.7 --useDCB --skipSystematics"
#     sed -i "s#'procs':'all',#'procs':'auto',#g" config_${year}_fiducial.py
# done
# python RunPackager.py --cats ggH0,ggH1,ggH2,ggH3,VBF0,VBF1,VBF2,VBF3 --exts fiducial_2016preVFP,fiducial_2016postVFP,fiducial_2017,fiducial_2018,fiducial_2022preEE,fiducial_2022postEE,fiducial_2023preBPix,fiducial_2023postBPix --mergeYears --batch local --outputExt packaged --outputPath ${outputpath}_proc
# python RunPlotter.py --procs all --cats all --years 2016preVFP,2016postVFP,2017,2018,2022preEE,2022postEE,2023preBPix,2023postBPix --ext packaged --outputPath ${outputpath}_proc
# for cat in ggH0 ggH1 ggH2 ggH3 VBF0 VBF1 VBF2 VBF3; do
#     echo "Plotting signal models for $cat category"
#     python RunPlotter.py --procs all --cats $cat --years 2016preVFP,2016postVFP,2017,2018,2022preEE,2022postEE,2023preBPix,2023postBPix --ext packaged --outputPath ${outputpath}_proc
# done

### Merge proc&year
# for cat in ggH0 ggH1 ggH2 ggH3 VBF0 VBF1 VBF2 VBF3; do
#     python scripts/calcPhotonSyst_sync.py --xvar mllg_cat_ --cat ${cat} --flavs ele,mu --ext ${outputExt} --scales CMS_scale_e,CMS_scale_g,CMS_scale_m --smears CMS_res_e,CMS_res_g,CMS_res_m
# done
# python RunSignalScripts.py --inputConfig config_flav_fiducial_nosys.py --mode signalFit --groupSignalFitJobsByCat --outputPath ${outputpath} --modeOpts "--ext ${outputExt} --doPlots --skipBeamspotReweigh --useDCB --skipSystematics" #  --skipSystematics
# python RunPackager.py --cats ggH0,ggH1,ggH2,ggH3,VBF0,VBF1,VBF2,VBF3 --exts ${outputExt} --mergeYears --batch local --outputExt packaged --outputPath ${outputpath} --massPoint 120,125,130
# python RunPlotter.py --procs all --cats ggH0,ggH1,ggH2,ggH3,VBF0,VBF1,VBF2,VBF3 --years all --ext packaged --outputPath ${outputpath}
# for cat in ggH0 ggH1 ggH2 ggH3 VBF0 VBF1 VBF2 VBF3; do
#     echo "Plotting signal models for $cat category"
#     python RunPlotter.py --procs all --cats $cat --years all --ext packaged --outputPath ${outputpath}
# done

# ggH0 ggH1 ggH2 ggH3 VBF0 VBF1 VBF2 VBF3 Incl0 Incl1 Incl2 Incl3 Incl4 Incl5 Incl6 Incl7
for cat in ggH0 ggH1 ggH2 ggH3 VBF0 VBF1 VBF2 VBF3; do 
    echo "generate signal combine pdf in $cat"
    python scripts/combineSignalPdf.py --cat $cat --outputPath ${outputpath} --inputExt packaged --outputExt combinedPDFs; 
done
