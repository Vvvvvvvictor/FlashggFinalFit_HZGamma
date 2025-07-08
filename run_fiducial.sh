#!/bin/bash

# Function to wait for Condor jobs to complete
wait_for_condor_jobs() {
    echo "Waiting for Condor jobs to complete..."
    local user=$(whoami)
    while true; do
        # Check for running, idle, or held jobs
        local running_jobs=$(condor_q $user -af JobStatus | grep -c '^[12]') # 1=Idle, 2=Running
        local held_jobs=$(condor_q $user -af JobStatus | grep -c '^5') # 5=Held

        if [ $held_jobs -gt 0 ]; then
            echo "$held_jobs jobs are held. Releasing them..."
            condor_release $user
            sleep 10 # Give some time for release command to take effect
            continue # Re-check immediately after releasing
        fi

        if [ $running_jobs -eq 0 ] && [ $held_jobs -eq 0 ]; then
            echo "All Condor jobs finished."
            break
        else
            local total_jobs=$((running_jobs + held_jobs))
            echo "Waiting for $total_jobs Condor jobs to finish... ($running_jobs running/idle, $held_jobs held). Checking again in 60 seconds."
            sleep 60
        fi
    done
}

clear
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsenv
echo $(pwd)
source $(pwd)/setup.sh
mainPath=$(pwd)
inputpath="input_final"
outputpath="output_final"

BackgroundWSPath="${mainPath}/Inputdata/${inputpath}/background"
BackgroundNtuplePath="${mainPath}/Inputdata/${inputpath}/fitting_bkg"
SignalNtuplePath="${mainPath}/Inputdata/${inputpath}/fitting_signal"
SignalWSPath="${mainPath}/Inputdata/${inputpath}/signal"
# BackgroundWSPath="${mainPath}/Inputdata/outputs_rzou_run3_finalfit/Background_workspace"
# BackgroundNtuplePath="${mainPath}/Inputdata/outputs_rzou_run3_finalfit/fitting_bkg"
# SignalNtuplePath="${mainPath}/Inputdata/outputs_rzou_run3_finalfit/fitting_signal"
# SignalWSPath="${mainPath}/Inputdata/outputs_rzou_run3_finalfit/Signal_workspace"

# SignalProcs=("ggH" "VBF" "WH" "ZH" "ttH") # "ggH" "VBF" "WplusH" "WminusH" "ZH" "ttH"
SignalProcs=("ggH" "VBF") # "ggH" "VBF"
mass_points=("125")
years=("2016preVFP" "2016postVFP" "2017" "2018" "2022preEE" "2022postEE" "2023preBPix" "2023postBPix") #"2016preVFP" "2016postVFP" "2017" "2018" "2022preEE" "2022postEE" "2023preBPix" "2023postBPix"
# years=("2017") #"2016preVFP" "2016postVFP" "2017" "2018" "2022preEE" "2022postEE" "2023preBPix" "2023postBPix"

############################################
# Tree2WS
############################################
cd ${mainPath}/Trees2WS
# make background ws
# mkdir -p ${BackgroundWSPath}
# python trees2ws_data.py --inputConfig config_Run2.py --inputTreeFile ${BackgroundNtuplePath}/Data/output_Data_all.root --outputWSDir ${BackgroundWSPath}

# # make signal ws
# mkdir -p ${SignalWSPath}
# for year in "${years[@]}"; do
#     for mass_point in "${mass_points[@]}"; do
#         for proc in "${SignalProcs[@]}"; do
#             mkdir -p ${SignalWSPath}_${year}
#             python trees2ws.py --inputConfig config_Run2.py --inputTreeFile ${SignalNtuplePath}/${proc}_M${mass_point}_${year}/output_${proc}_M${mass_point}.root --inputMass ${mass_point} --productionMode ${proc} --year ${year} --outputWSDir ${SignalWSPath}_${year} #--doSystematics
#         done
#     done
# done

# for year in "${years[@]}"; do
#     for proc in "${SignalProcs[@]}"; do
#         python mass_shifter.py --inputMass 125 --targetMass 120 --inputWSFile ${SignalWSPath}_${year}/output_${proc}_M125_pythia8_${proc}.root
#         python mass_shifter.py --inputMass 125 --targetMass 130 --inputWSFile ${SignalWSPath}_${year}/output_${proc}_M125_pythia8_${proc}.root
#     done
# done

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

# signalfit
for year in "${years[@]}"; do
    echo "Running signal fitting for year $year"
    # sed -i "s#input_jiehan#input_final#g" config_${year}_fiducial.py
    python RunSignalScripts.py --inputConfig config_${year}_fiducial.py --mode signalFit --groupSignalFitJobsByCat --outputPath ${outputpath} --modeOpts "--doPlots --beamspotWidthData 3.5 --beamspotWidthMC 3.7 --useDCB --skipSystematics"
done

# # packaged
# echo "Packaging signal models"
python RunPackager.py --cats ggH0,ggH1,ggH2,ggH3,VBF0,VBF1,VBF2,VBF3 --exts fiducial_2016preVFP,fiducial_2016postVFP,fiducial_2017,fiducial_2018,fiducial_2022preEE,fiducial_2022postEE,fiducial_2023preBPix,fiducial_2023postBPix --mergeYears --batch local --outputExt packaged --outputPath ${outputpath}
# ggH0,ggH1,ggH2,ggH3,VBF0,VBF1,VBF2,VBF3,VHlep,ZHinv,ttHl,ttHh
# Incl0,Incl1,Incl2,Incl3,Incl4,Incl5,Incl6,Incl7

# # signal model plotting
echo "Plotting signal models"
python RunPlotter.py --procs all --cats all --years 2016preVFP,2016postVFP,2017,2018,2022preEE,2022postEE,2023preBPix,2023postBPix --ext packaged --outputPath ${outputpath}

for proc in ggH VBF; do
    echo "Plotting signal models for $proc production mode"
    python RunPlotter.py --procs $proc --cats all --years 2016preVFP,2016postVFP,2017,2018,2022preEE,2022postEE,2023preBPix,2023postBPix --ext packaged --outputPath ${outputpath}
done

for cat in ggH0 ggH1 ggH2 ggH3 VBF0 VBF1 VBF2 VBF3; do
    echo "Plotting signal models for $cat category"
    python RunPlotter.py --procs all --cats $cat --years 2016preVFP,2016postVFP,2017,2018,2022preEE,2022postEE,2023preBPix,2023postBPix --ext packaged --outputPath ${outputpath}
done

for flav in ele mu; do
    echo "Plotting signal models for $flav lepton flavor"
    python RunPlotter.py --procs all --cats all --flavs $flav --years 2016preVFP,2016postVFP,2017,2018,2022preEE,2022postEE,2023preBPix,2023postBPix --ext packaged --outputPath ${outputpath}
done

# python RunPlotter.py --procs all --cats ggH0,ggH1,ggH2,ggH3 --years 2016preVFP,2016postVFP,2017,2018,2022preEE,2022postEE,2023preBPix,2023postBPix --ext packaged
# python RunPlotter.py --procs all --cats all --years 2016preVFP --ext packaged

# ggH0 ggH1 ggH2 ggH3 VBF0 VBF1 VBF2 VBF3 Incl0 Incl1 Incl2 Incl3 Incl4 Incl5 Incl6 Incl7
# for cat in ggH0 ggH1 ggH2 ggH3; do python scripts/combineSignalPdf.py --cat $cat; done
for cat in ggH0 ggH1 ggH2 ggH3 VBF0 VBF1 VBF2 VBF3; do 
    echo "generate signal combine pdf in $cat"
    python scripts/combineSignalPdf.py --cat $cat --outputPath ${outputpath}; 
done

###########################################
# Background
###########################################
cd ${mainPath}/Background
# make clean; make -j 16;
# rm -rf outdir_rzou_run3/* dat/*
# rm -rf outdir_jiehan/* dat/*
# rm -rf params_rzou_run3/*
# rm -rf params_jiehan/*
# python RunBackgroundScripts.py --inputConfig config_fiducial_run2.py --mode fTestParallel #--jobOpts “--blindFit”
# mv params/ params_rzou_run3/
# mv params/ params_jiehan/

## Sync
# make workspace with documentation
# ./bin/BkgWSMaker -t "pow1,exp3,lau2" -i /eos/user/m/mingtao/workspace/zgamma/CMSSW_10_2_13/src/FlashggFinalFit_HZGamma/Inputdata/outputs_rzou_run3_finalfit/Background_workspace/ws/output_Data_all.root --saveMultiPdf outdir_fiducialAnalysis_sync/CMS-HGG_multipdf_ggH0.root -D outdir_fiducialAnalysis_sync/bkgfTest-Data -f ggH0 --mgg_low 105  --isData 1 --year all --catOffset 0 -v
# ./bin/BkgWSMaker -t "bern4,pow3,exp3,lau3" -i /eos/user/m/mingtao/workspace/zgamma/CMSSW_10_2_13/src/FlashggFinalFit_HZGamma/Inputdata/outputs_rzou_run3_finalfit/Background_workspace/ws/output_Data_all.root --saveMultiPdf outdir_fiducialAnalysis_sync/CMS-HGG_multipdf_ggH1.root -D outdir_fiducialAnalysis_sync/bkgfTest-Data -f ggH1 --mgg_low 105  --isData 1 --year all --catOffset 1 -v
# ./bin/BkgWSMaker -t "bern5,pow3,exp3,lau4" -i /eos/user/m/mingtao/workspace/zgamma/CMSSW_10_2_13/src/FlashggFinalFit_HZGamma/Inputdata/outputs_rzou_run3_finalfit/Background_workspace/ws/output_Data_all.root --saveMultiPdf outdir_fiducialAnalysis_sync/CMS-HGG_multipdf_ggH2.root -D outdir_fiducialAnalysis_sync/bkgfTest-Data -f ggH2 --mgg_low 100  --isData 1 --year all --catOffset 2 -v
# ./bin/BkgWSMaker -t "bern4,pow3,exp3,lau2,modgau1" -i /eos/user/m/mingtao/workspace/zgamma/CMSSW_10_2_13/src/FlashggFinalFit_HZGamma/Inputdata/outputs_rzou_run3_finalfit/Background_workspace/ws/output_Data_all.root --saveMultiPdf outdir_fiducialAnalysis_sync/CMS-HGG_multipdf_ggH3.root -D outdir_fiducialAnalysis_sync/bkgfTest-Data -f ggH3 --mgg_low 95  --isData 1 --year all --catOffset 3 -v


###########################################
# Datacard
###########################################
cd ${mainPath}/Datacard
# for cat in "ggH0" "ggH1" "ggH2" "ggH3" "VBF0" "VBF1" "VBF2" "VBF3"; do # "ggH0" "ggH1" "ggH2" "ggH3" "VBF0" "VBF1" "VBF2" "VBF3"
#     inputstring=$(IFS=,; for year in "${years[@]}"; do echo -n "${year}=${SignalWSPath}_${year}/,"; done | sed 's/,$//')
#     echo $inputstring
#     python RunYields.py --inputWSDirMap $inputstring --cats $cat --procs auto --batch local --ext $cat --mergeYears --skipZeroes # --doSystematics
#     python makeDatacard.py --years $(IFS=,; echo "${years[*]}") --ext $cat --prune --pruneThreshold 0.001 --mergeYears --output Datacard_$cat #--doSystematics 
# done

# inputstring=$(IFS=,; for year in "${years[@]}"; do echo -n "${year}=${SignalWSPath}_${year}/,"; done | sed 's/,$//')
# python RunYields.py --inputWSDirMap $inputstring --cats auto --procs auto --batch local --ext fiducial --mergeYears --skipZeroes # --doSystematics
# python makeDatacard.py --years $(IFS=,; echo "${years[*]}") --ext fiducial --prune --pruneThreshold 0.001 --mergeYears --output Datacard_fiducial #--doSystematics

# python RunYields.py --inputWSDirMap $inputstring --cats ggH0,ggH1,ggH2,ggH3 --procs auto --doSystematics --batch local --ext ggH --mergeYears --skipZeroes
# python makeDatacard.py --years $(IFS=,; echo "${years[*]}") --ext ggH --prune --pruneThreshold 0.001 --doSystematics --mergeYears --output Datacard_ggH
# python RunYields.py --inputWSDirMap $inputstring --cats VBF0,VBF1,VBF2,VBF3 --procs auto --doSystematics --batch local --ext VBF --mergeYears --skipZeroes
# python makeDatacard.py --years $(IFS=,; echo "${years[*]}") --ext VBF --prune --pruneThreshold 0.001 --doSystematics --mergeYears --output Datacard_VBF
# python RunYields.py --inputWSDirMap $inputstring --cats VHlep,ZHinv,ttHl,ttHh --procs auto --doSystematics --batch local --ext others --mergeYears --skipZeroes
# python makeDatacard.py --years $(IFS=,; echo "${years[*]}") --ext others --prune --pruneThreshold 0.001 --doSystematics --mergeYears --output Datacard_others

###########################################
# Combine
###########################################
echo Combine
cd ${mainPath}/Combine
# cp ../Datacard/Datacard*.txt .
# rm -rf Models; mkdir Models
# mkdir Models/signal
# mkdir Models/background
# cp -r ../Signal/${outputpath}/outdir_packaged/CMS-HGG_sigfit_packaged_* Models/signal/
# cp -r ../Background/outdir_final/CMS-HGG_multipdf_* Models/background/

# for ext in "_ggH0" "_ggH1" "_ggH2" "_ggH3" "_VBF0" "_VBF1"  "_VBF2" "_VBF3"; do # "_fiducial" "_ggH0" "_ggH1" "_ggH2" "_ggH3" "_VBF0" "_VBF1"  "_VBF2" "_VBF3"; do
#     # python RunText2Workspace.py --ext $ext --mode mu_fiducial --batch local
#     # combine Datacard${ext}_mu_fiducial.root -M Significance -t -1 --expectSignal=1 -m 125.0 -n $ext --setParameters MH=125,pdfindex_VBF1_13TeV=2,pdfindex_VBF2_13TeV=1 --verbose 9 >> log${ext}
#     combineTool.py -d Datacard${ext}_mu_fiducial.root -n ${ext} -M MultiDimFit --algo grid --points 50 --noMCbonly 1  --cminDefaultMinimizerStrategy 0 --cminApproxPreFitTolerance=100  --cminFallbackAlgo Minuit2,Migrad,0:0.1 --cminDefaultMinimizerTolerance 0.1 --X-rtd MINIMIZER_MaxCalls=9999999  --X-rtd MINIMIZER_analytic --X-rtd FAST_VERTICAL_MORPH  --X-rtd MINIMIZER_freezeDisassociatedParams  --X-rtd OPTIMIZE_BOUNDS=0  -v 3  --setParameterRanges r=-1,3 -m 125  --floatOtherPOIs 1 --X-rtd NO_INITIAL_SNAP  --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2  --saveWorkspace --alignEdges 1 --toysFrequentist -t -1 --expectSignal=1 >> log_scan${ext}
# done

# # Bias Study
# cd Checks/
# # Only to make directories
# for cat in "Incl0"; do # "VBF0" "VBF1" "VBF2" "VBF3"
#     python RunBiasStudy.py -j $cat -d ../Datacard_${cat}_mu_fiducial.root -t -n 1000 --dryRun
#     python RunBiasStudy.py -j $cat -d ../Datacard_${cat}_mu_fiducial.root -f -n 1000 -c "--cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints   --rMin -100 --rMax 100  --freezeParameters MH" --dryRun
#     python RunBiasStudy.py -j $cat -d ../Datacard_${cat}_mu_fiducial.root -f -n 100 -c "--alignEdges 1 --setParameterRanges CMS_hgg_mass=110,130 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2"
#     python RunBiasStudy.py -j $cat -d ../Datacard_${cat}_mu_fiducial.root -p --gaussianFit
# done

# python Checks/RunBiasStudy.py -d /afs/cern.ch/user/j/jiehan/finalfit/CMSSW_10_2_13/src/flashggFinalFit/Combine/Datacard_VBF0_mu_fiducial.root -t 
# python Checks/RunBiasStudy.py -d /afs/cern.ch/user/j/jiehan/finalfit/CMSSW_10_2_13/src/flashggFinalFit/Combine/Datacard_VBF0_mu_fiducial.root -f -c "--cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints   --rMin -100 --rMax 100  --freezeParameters MH"
# python Checks/RunBiasStudy.py -d /afs/cern.ch/user/j/jiehan/finalfit/CMSSW_10_2_13/src/flashggFinalFit/Combine/Datacard_VBF0_mu_fiducial.root -p --gaussianFit

# python Checks/RunBiasStudy.py -d /afs/cern.ch/user/j/jiehan/finalfit/CMSSW_10_2_13/src/flashggFinalFit/Combine/Datacard_fiducial_mu_fiducial.root -t -f -p --gaussianFit -c "--cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints   --rMin -100 --rMax 100  --freezeParameters MH"

# 改改名字
# sed -i -E '/13TeV_bkgshape/s/(13TeV_bkgshape)/2022_\1_norm/g; /pdfindex/s/(resolution)/\1_2022/g' Datacard.txt
# # 送到afs去
# cp models.py ~/finalfit/CMSSW_10_2_13/src/flashggFinalFit/Combine/
# cp inputs.json ~/finalfit/CMSSW_10_2_13/src/flashggFinalFit/Combine/
# cp Datacard*.txt ~/finalfit/CMSSW_10_2_13/src/flashggFinalFit/Combine/
# cp -r Models ~/finalfit/CMSSW_10_2_13/src/flashggFinalFit/Combine/
# cp Checks/RunBiasStudy.py ~/finalfit/CMSSW_10_2_13/src/flashggFinalFit/Combine/Checks/
# cp Checks/biasUtils.py ~/finalfit/CMSSW_10_2_13/src/flashggFinalFit/Combine/Checks/
# cd ~/finalfit/CMSSW_10_2_13/src/flashggFinalFit/
# # 环境
# source /cvmfs/cms.cern.ch/cmsset_default.sh
# cmsenv
# source setup.sh
# cd Combine
# python RunText2Workspace.py --ext _fiducial --mode mu_fiducial --batch local
# # python RunFits.py --inputJson inputs.json --ext _fiducial --mode mu_fiducial --batch condor

# for cat in "VBF0" "VBF1" "VBF2" "VBF3"; do
#     python RunText2Workspace.py --ext _$cat --mode mu_fiducial --batch local
#     # python RunFits.py --inputJson inputs.json --ext _$cat --mode mu_fiducial --batch condor
# done

# # # for ext in "_VBF"; #"_fiducial" "_ggH" "_VBF" "_others"; 
# # # do
# # #     python RunText2Workspace.py --ext $ext --mode mu_fiducial --batch local
# # #     # python RunFits.py --inputJson inputs.json --ext $ext --mode mu_fiducial --batch condor
# # # done

# # # for ext in "_fiducial" "_VBF0" "_VBF1" "_VBF2" "_VBF3"; do cd runFits${ext}_mu_fiducial; rm *log *err *out *root; find . -name "*statonly*.sub" -exec condor_submit {} \;; cd ..; done # "_ggH" "_VBF" "_others"

# # Bias Study
# cd Checks/
# for cat in "VBF0" "VBF1"  "VBF2" "VBF3"; do # 
#     python RunBiasStudy.py -j $cat -d /eos/home-j/jiehan/finalfit_102X/CMSSW_10_2_13/src/flashggFinalFit/Combine/Datacard_${cat}_mu_fiducial.root -t -n 1000 --condor --dryRun
#     python RunBiasStudy.py -j $cat -d /eos/home-j/jiehan/finalfit_102X/CMSSW_10_2_13/src/flashggFinalFit/Combine/Datacard_${cat}_mu_fiducial.root -f -n 1000 -c "--cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints   --rMin -100 --rMax 100  --freezeParameters MH" --condor --dryRun
# done
# # need to submit with condor and better to sleep 600s between two jobs to avoid overflow of afs quota
# find . -maxdepth 1 -name '*toys.sub' -exec condor_submit {} \; -exec echo 'sleep 10 minutes...' \; -exec sleep 600 \;
# find . -maxdepth 1 -name '*fits.sub' -exec condor_submit {} \; -exec echo 'sleep 10 minutes...' \; -exec sleep 600 \;

# # Call the function to wait for jobs
# wait_for_condor_jobs

# # Now proceed with the next steps after jobs are done
# echo "Proceeding with plotting bias results..."
# for cat in "VBF0" "VBF1" "VBF2" "VBF3"; do
#     python RunBiasStudy.py -j $cat -d /eos/home-j/jiehan/finalfit_102X/CMSSW_10_2_13/src/flashggFinalFit/Combine/Datacard_${cat}_mu_fiducial.root -p --gaussianFit
# done

# cd .. # Go back to Combine directory from Checks

# for cat in "fiducial" "VBF0" "VBF1" "VBF2" "VBF3"; do
#     combine Datacard_${cat}_mu_fiducial.root -M Significance -t -1 --expectSignal=1 -m 125.0 -n allCats
# done

# # submit mu scan condor jobs
# for ext in "_fiducial" "_VBF0" "_VBF1" "_VBF2" "_VBF3"; do cd runFits${ext}_mu_fiducial; find . -name "*.sub" -exec condor_submit {} \;; cd ..; done

# # Call the function to wait for jobs
# wait_for_condor_jobs

# python CollectFits.py --inputJson inputs.json --mode mu_fiducial --ext _fiducial
# plot1DScan.py runFits_fiducial_mu_fiducial/profile1D_syst_r.root --y-cut 5 --y-max 5 --output r_fiducial_fixed_statsyst --POI r --translate ../Plots/pois_mu.json --main-label 'Expected' --main-color 1 --others runFits_fiducial_mu_fiducial/profile1D_statonly_r.root:'Stat only':2 --logo-sub 'Work in Progress'

# # exts=("_fiducial" "_ggH" "_VBF" "_others")
# # lims=("5" "2" "3" "0.2")
# # for i in {0..3}; do
# #     ext=${exts[i]}
# #     lim=${lims[i]}
# #     python CollectFits.py --inputJson inputs.json --mode mu_fiducial --ext $ext
# #     plot1DScan.py runFits${ext}_mu_fiducial/profile1D_syst_r.root --y-cut $lim --y-max $lim --output r${ext}_fixed_statsyst --POI r --translate ../Plots/pois_mu.json --main-label 'Expected' --main-color 1 --others runFits${ext}_mu_fiducial/profile1D_statonly_r.root:'Stat only':2 --logo-sub 'Work in Progress'
# # done

# for cat in "VBF0" "VBF1" "VBF2" "VBF3"; do
#     python CollectFits.py --inputJson inputs.json --mode mu_fiducial --ext _$cat
#     plot1DScan.py runFits_${cat}_mu_fiducial/profile1D_syst_r.root --y-cut 1 --y-max 1 --output r_${cat}_fiducial_fixed_statsyst --POI r --translate ../Plots/pois_mu.json --main-label 'Expected' --main-color 1 --others runFits_${cat}_mu_fiducial/profile1D_statonly_r.root:'Stat only':2 --logo-sub 'Work in Progress'
# done

# cd ../Plots
# python makeSplusBModelPlot.py --inputWSFile ../Combine/Datacard_fiducial_mu_fiducial.root --cats all --doZeroes --ext _Run23 --translateCats cats.json

# impact plot
# text2workspace.py Datacard_fiducial.txt -m 125 -o Datacard_fiducial.root
# combine -M AsymptoticLimits -m 125 -n MX3000_MH125 -d MX3000_MH125.txt --run expected --freezeParameters MH
# combineTool.py -M Impacts -d Datacard_fiducial.root -m 125 --freezeParameters MH -n .impacts --setParameterRanges r=0,2 --doInitialFit -t -1 --expectSignal 1 --robustFit 1 --X-rt MINIMIZER_freezeDisassociatedParams   --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --cminDefaultMinimizerStrategy 0 
# combineTool.py -M Impacts -d Datacard_fiducial.root -m 125 --freezeParameters MH -n .impacts --setParameterRanges r=0,2 --doFits --robustFit 1 -t -1 --expectSignal 1   --X-rt MINIMIZER_freezeDisassociatedParams   --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --cminDefaultMinimizerStrategy 0 
# combineTool.py -M Impacts -d Datacard_fiducial.root -m 125 --freezeParameters MH -n .impacts --setParameterRanges r=0,2 -o Datacard_fiducial.json -t -1 --X-rt MINIMIZER_freezeDisassociatedParams   --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --cminDefaultMinimizerStrategy 0 
# plotImpacts.py -i Datacard_fiducial.json -o Datacard_fiducial

# 定义要执行的命令列表
commands=(
    # # Signal
    # "cd ${mainPath}/Signal"
    # # fTest
    # "python RunSignalScripts.py --inputConfig config_2022preEE_fiducial.py --mode fTest --modeOpts "--doPlots""
    # "python RunSignalScripts.py --inputConfig config_2022postEE_fiducial.py --mode fTest --modeOpts "--doPlots""
    # # syst
    # "python RunSignalScripts.py --inputConfig config_2022preEE_fiducial.py --mode calcPhotonSyst"
    # "python RunSignalScripts.py --inputConfig config_2022postEE_fiducial.py --mode calcPhotonSyst"
    # # signalfit
    # # "python RunSignalScripts.py --inputConfig config_2022preEE_fiducial.py --mode signalFit --groupSignalFitJobsByCat --modeOpts "--doPlots --beamspotWidthData 3.5 --beamspotWidthMC 3.7" "
    # "python RunSignalScripts.py --inputConfig config_2022preEE_fiducial.py --mode signalFit --groupSignalFitJobsByCat --modeOpts "--doPlots --beamspotWidthData 3.5 --beamspotWidthMC 3.7" "
    # "python RunSignalScripts.py --inputConfig config_2022postEE_fiducial.py --mode signalFit --groupSignalFitJobsByCat --modeOpts "--doPlots --beamspotWidthData 3.5 --beamspotWidthMC 3.7" "
    # # packaged
    # "python RunPackager.py --cats best_resolution,medium_resolution,worst_resolution --exts fiducial_2022preEE,fiducial_2022postEE --mergeYears --batch local"
    # Background
    # "cd ${mainPath}/Background"
    # "python RunBackgroundScripts.py --inputConfig config_2022_fiducial.py --mode fTestParallel"
    # # Datacard
    # "cd ../Datacard"
    # "python RunYields.py --inputWSDirMap 2022preEE=/eos/user/c/chpan/input_finalfit/signal_fiducial_2022preEE_sublast/,2022postEE=/eos/user/c/chpan/input_finalfit/signal_fiducial_2022postEE_sublast/ --cats auto --procs auto --doSystematics --batch local --ext fiducial --mergeYears --skipZeroes"
    # "python makeDatacard.py --years 2022preEE,2022postEE --ext fiducial --prune --pruneThreshold 0.001 --doSystematics --doMCStatUncertainty"
    # # Combine
    # "cd ${mainPath}/Combine"
    # "cp ../Datacard/Datacard_fiducial.txt ."
    # "rm -rf Models; mkdir Models"
    # "mkdir Models/signal"
    # "mkdir Models/background"
    # "cp -r ../Signal/outdir_packaged/CMS-HGG_sigfit_packaged_* Models/signal/"
    # "cp -r ../Background/outdir_fiducialAnalysis/CMS-HGG_multipdf_* Models/background/"
    # # # 改改名字
    # # "sed -i -E '/13TeV_bkgshape/s/(13TeV_bkgshape)/2022_\1_norm/g; /pdfindex/s/(resolution)/\1_2022/g' Datacard.txt"
    # # 送到afs去
    # "cp input.json ~/finalfit/CMSSW_10_2_13/src/flashggFinalFit/Combine/"
    # "cp Datacard_fiducial.txt ~/finalfit/CMSSW_10_2_13/src/flashggFinalFit/Combine/"
    # "rm -rf ~/finalfit/CMSSW_10_2_13/src/flashggFinalFit/Combine/Models"
    # "cp -r Models ~/finalfit/CMSSW_10_2_13/src/flashggFinalFit/Combine/"
    # "cd ~/finalfit/CMSSW_10_2_13/src/flashggFinalFit/"
    # 环境
    # "source /cvmfs/cms.cern.ch/cmsset_default.sh"
    # "cmsenv"
    # "source setup.sh"
    # Combineinputs.json
    # "cd Combine"
    # "python RunText2Workspace.py --ext _fiducial --mode mu_fiducial --batch local"
    # "python RunFits.py --inputJson inputs.json --mode mu_fiducial --batch condor"
    # # 在这里添加更多命令
)

# 遍历命令列表，依次执行每个命令
# for cmd in "${commands[@]}"; do
#     echo "Executing command: $cmd"
#     # 执行命令
#     eval "$cmd"
#     # 检查命令执行状态
#     if [ $? -eq 0 ]; then
#         echo "Command executed successfully"
#     else
#         echo "Error executing command: $cmd"
#         # 可选择在错误发生时终止执行
#         # exit 1
#     fi
# done

# # 在此处添加代码以等待 Condor 作业完成#     eval "$cmd"n#     # 检查命令执行状态h#     if [ $? -eq 0 ]; then #         echo "Command executed successfully"e#     else #         echo "Error executing command: $cmd" #         # 可选择在错误发生时终止执行 #         # exit 1s#     fis# doneo echo "All commands executed"