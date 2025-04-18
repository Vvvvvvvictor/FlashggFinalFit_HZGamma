#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsenv
source setup.sh
mainPath=$(pwd)
BackgroundWSPath="/eos/home-j/jiehan/root/input_finalfit/background"
SignalNtuplePath="/eos/home-j/jiehan/root/fitting_signal"
SignalWSPath="/eos/home-j/jiehan/root/input_finalfit/signal"

SinalProcs=("ggH" "VBF" "WH" "ZH" "ttH") # "ggH" "VBF" "WplusH" "WminusH" "ZH" "ttH"
mass_points=("125")
years=("2022preEE" "2022postEE" "2023preBPix" "2023postBPix")

############################################
# Tree2WS
############################################
cd ${mainPath}/Trees2WS
# make background ws
mkdir -p ${BackgroundWSPath}
python trees2ws_data.py --inputConfig config_Run3.py --inputTreeFile /eos/home-j/jiehan/root/fitting_bkg/Data/output_Data_Run3.root --outputWSDir ${BackgroundWSPath}

# make signal ws
mkdir -p ${SignalWSPath}
for year in "${years[@]}"; do
    for mass_point in "${mass_points[@]}"; do
        for proc in "${SinalProcs[@]}"; do
            mkdir -p ${SignalWSPath}_${year}
            python trees2ws.py --inputConfig config_Run3.py --inputTreeFile ${SignalNtuplePath}/${proc}_M${mass_point}_${year}/output_${proc}_M${mass_point}.root --inputMass ${mass_point} --productionMode ${proc} --year ${year} --outputWSDir ${SignalWSPath}_${year} #--doSystematics
        done
    done
done

for year in "${years[@]}"; do
    for proc in "${SinalProcs[@]}"; do
        python mass_shifter.py --inputMass 125 --targetMass 120 --inputWSFile ${SignalWSPath}_${year}/output_${proc}_M125_pythia8_${proc}.root
        python mass_shifter.py --inputMass 125 --targetMass 130 --inputWSFile ${SignalWSPath}_${year}/output_${proc}_M125_pythia8_${proc}.root
    done
done

#########################################
# Signal
#########################################
cd ${mainPath}/Signal
# # fTest
# for year in "${years[@]}"; do
#     python RunSignalScripts.py --inputConfig config_${year}_fiducial.py --mode fTest --modeOpts "--doPlots"
# done
# # syst
# for year in "${years[@]}"; do
#     python RunSignalScripts.py --inputConfig config_${year}_fiducial.py --mode calcPhotonSyst
# done
# signalfit
for year in "${years[@]}"; do
    python RunSignalScripts.py --inputConfig config_${year}_fiducial.py --mode signalFit --groupSignalFitJobsByCat --modeOpts "--doPlots --beamspotWidthData 3.5 --beamspotWidthMC 3.7 --useDCB --skipSystematics" 
done
# signal model plotting
python RunPlotter.py --procs all --cats all --years 2022preEE,2022postEE,2023preBPix,2023postBPix --ext packaged

# packaged
python RunPackager.py --cats VBF0,VBF1,VBF2,VBF3 --exts fiducial_2022preEE,fiducial_2022postEE,fiducial_2023preBPix,fiducial_2023postBPix --mergeYears --batch local --outputExt packaged_run3

# for cat in VBF0 VBF1 VBF2 VBF3; do python scripts/combineSignalPdf.py --cat $cat; done

###########################################
# Background
###########################################
cd ${mainPath}/Background
# make clean; make -j 16;
python RunBackgroundScripts.py --inputConfig config_fiducial_run3.py --mode fTestParallel --jobOpts “--blindFit”

###########################################
# Datacard
###########################################
cd ${mainPath}/Datacard
for cat in "VBF0" "VBF1" "VBF2" "VBF3"; do
    inputstring=$(IFS=,; for year in "${years[@]}"; do echo -n "${year}=${SignalWSPath}_${year}/,"; done | sed 's/,$//')
    echo $inputstring
    python RunYields.py --inputWSDirMap $inputstring --cats $cat --procs auto --batch local --ext ${cat}_run3 --mergeYears --skipZeroes # --doSystematics
    python makeDatacard.py --years $(IFS=,; echo "${years[*]}") --ext ${cat}_run3 --prune --pruneThreshold 0.001 --mergeYears --output Datacard_${cat}_run3 #--doSystematics 
done
inputstring=$(IFS=,; for year in "${years[@]}"; do echo -n "${year}=${SignalWSPath}_${year}/,"; done | sed 's/,$//')
python RunYields.py --inputWSDirMap $inputstring --cats auto --procs auto --batch local --ext fiducial_run3 --mergeYears --skipZeroes # --doSystematics
python makeDatacard.py --years $(IFS=,; echo "${years[*]}") --ext fiducial_run3 --prune --pruneThreshold 0.001 --mergeYears --output Datacard_fiducial_run3 #--doSystematics

# # python RunYields.py --inputWSDirMap $inputstring --cats ggH0,ggH1,ggH2,ggH3 --procs auto --doSystematics --batch local --ext ggH --mergeYears --skipZeroes
# # python makeDatacard.py --years $(IFS=,; echo "${years[*]}") --ext ggH --prune --pruneThreshold 0.001 --doSystematics --mergeYears --output Datacard_ggH
# # python RunYields.py --inputWSDirMap $inputstring --cats VBF0,VBF1,VBF2,VBF3 --procs auto --doSystematics --batch local --ext VBF --mergeYears --skipZeroes
# # python makeDatacard.py --years $(IFS=,; echo "${years[*]}") --ext VBF --prune --pruneThreshold 0.001 --doSystematics --mergeYears --output Datacard_VBF
# # python RunYields.py --inputWSDirMap $inputstring --cats VHlep,ZHinv,ttHl,ttHh --procs auto --doSystematics --batch local --ext others --mergeYears --skipZeroes
# # python makeDatacard.py --years $(IFS=,; echo "${years[*]}") --ext others --prune --pruneThreshold 0.001 --doSystematics --mergeYears --output Datacard_others

###########################################
# Combine
###########################################
# cd ${mainPath}/Combine
# cp ../Datacard/Datacard*.txt .
# rm -rf Models; mkdir Models
# mkdir Models/signal
# mkdir Models/background
# cp -r ../Signal/outdir_packaged/CMS-HGG_sigfit_packaged_* Models/signal/
# cp -r ../Background/outdir_fiducialAnalysis/CMS-HGG_multipdf_* Models/background/

# for ext in "_fiducial" "_VBF0" "_VBF1"  "_VBF2" "_VBF3"; do #"_fiducial" "_VBF0" "_VBF1"  "_VBF2" "_VBF3"; do
#     python RunText2Workspace.py --ext $ext --mode mu_fiducial --batch local
# done

# # Bias Study
# cd Checks/
# for cat in "VBF0" "VBF1"  "VBF2" "VBF3"; do
#     # python RunBiasStudy.py -j $cat -d ../Datacard_${cat}_mu_fiducial.root -t -e 0 -n 1000
#     # python RunBiasStudy.py -j $cat -d ../Datacard_${cat}_mu_fiducial.root -f -e 0 -n 1000 -c"--cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints   --rMin -100 --rMax 100  --freezeParameters MH"
#     # python RunBiasStudy.py -j $cat -d ../Datacard_${cat}_mu_fiducial.root -f -n 100 -c "--alignEdges 1 --setParameterRanges CMS_hgg_mass=110,130 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2"
#     python RunBiasStudy.py -j $cat -e 0 -d ../Datacard_${cat}_mu_fiducial.root -p --gaussianFit
# done
# python Checks/RunBiasStudy.py -d /afs/cern.ch/user/j/jiehan/finalfit/CMSSW_10_2_13/src/flashggFinalFit/Combine/Datacard_VBF0_mu_fiducial.root -t 
# python Checks/RunBiasStudy.py -d /afs/cern.ch/user/j/jiehan/finalfit/CMSSW_10_2_13/src/flashggFinalFit/Combine/Datacard_VBF0_mu_fiducial.root -f -c "--cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints   --rMin -100 --rMax 100  --freezeParameters MH"
# python Checks/RunBiasStudy.py -d /afs/cern.ch/user/j/jiehan/finalfit/CMSSW_10_2_13/src/flashggFinalFit/Combine/Datacard_VBF0_mu_fiducial.root -p --gaussianFit

# python Checks/RunBiasStudy.py -d /afs/cern.ch/user/j/jiehan/finalfit/CMSSW_10_2_13/src/flashggFinalFit/Combine/Datacard_fiducial_mu_fiducial.root -t -f -p --gaussianFit -c "--cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints   --rMin -100 --rMax 100  --freezeParameters MH"

# 改改名字
# sed -i -E '/13TeV_bkgshape/s/(13TeV_bkgshape)/2022_\1_norm/g; /pdfindex/s/(resolution)/\1_2022/g' Datacard.txt
# 送到afs去
cp models.py ~/finalfit/CMSSW_10_2_13/src/flashggFinalFit/Combine/
cp inputs.json ~/finalfit/CMSSW_10_2_13/src/flashggFinalFit/Combine/
cp Datacard*.txt ~/finalfit/CMSSW_10_2_13/src/flashggFinalFit/Combine/
cp -r Models ~/finalfit/CMSSW_10_2_13/src/flashggFinalFit/Combine/
cd ~/finalfit/CMSSW_10_2_13/src/flashggFinalFit/
# 环境
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsenv
source setup.sh
cd Combine
python RunText2Workspace.py --ext _fiducial --mode mu_fiducial --batch local
# python RunFits.py --inputJson inputs.json --ext _fiducial --mode mu_fiducial --batch condor

# # # for ext in "_VBF"; #"_fiducial" "_ggH" "_VBF" "_others"; 
# # # do
# # #     python RunText2Workspace.py --ext $ext --mode mu_fiducial --batch local
# # #     # python RunFits.py --inputJson inputs.json --ext $ext --mode mu_fiducial --batch condor
# # # done

# # # for ext in "_fiducial" "_VBF0" "_VBF1" "_VBF2" "_VBF3"; do cd runFits${ext}_mu_fiducial; rm *log *err *out *root; find . -name "*statonly*.sub" -exec condor_submit {} \;; cd ..; done # "_ggH" "_VBF" "_others"

for cat in "VBF0" "VBF1" "VBF2" "VBF3"; do
    python RunText2Workspace.py --ext _$cat --mode mu_fiducial --batch local
    # python RunFits.py --inputJson inputs.json --ext _$cat --mode mu_fiducial --batch condor
done

# for ext in "_fiducial" "_VBF0" "_VBF1" "_VBF2" "_VBF3"; do cd runFits${ext}_mu_fiducial; find . -name "*.sub" -exec condor_submit {} \;; cd ..; done

# python CollectFits.py --inputJson inputs.json --mode mu_fiducial --ext _fiducial
# plot1DScan.py runFits_fiducial_mu_fiducial/profile1D_syst_r.root --y-cut 5 --y-max 5 --output r_fiducial_fixed_statsyst --POI r --translate ../Plots/pois_mu.json --main-label 'Expected' --main-color 1 --others runFits_fiducial_mu_fiducial/profile1D_statonly_r.root:'Stat only':2 --logo-sub 'Work in Progress'

# exts=("_fiducial" "_ggH" "_VBF" "_others")
# lims=("5" "2" "3" "0.2")
# for i in {0..3}; do
#     ext=${exts[i]}
#     lim=${lims[i]}
#     python CollectFits.py --inputJson inputs.json --mode mu_fiducial --ext $ext
#     plot1DScan.py runFits${ext}_mu_fiducial/profile1D_syst_r.root --y-cut $lim --y-max $lim --output r${ext}_fixed_statsyst --POI r --translate ../Plots/pois_mu.json --main-label 'Expected' --main-color 1 --others runFits${ext}_mu_fiducial/profile1D_statonly_r.root:'Stat only':2 --logo-sub 'Work in Progress'
# done

# for cat in "VBF0" "VBF1" "VBF2" "VBF3"; do
#     python CollectFits.py --inputJson inputs.json --mode mu_fiducial --ext _$cat
#     # plot1DScan.py runFits_${cat}_mu_fiducial/profile1D_syst_r.root --y-cut 1 --y-max 1 --output r_${cat}_fiducial_fixed_statsyst --POI r --translate ../Plots/pois_mu.json --main-label 'Expected' --main-color 1 --others runFits_${cat}_mu_fiducial/profile1D_statonly_r.root:'Stat only':2 --logo-sub 'Work in Progress'
# done

# cd ../Plots
# python makeSplusBModelPlot.py --inputWSFile ../Combine/Datacard_fiducial_mu_fiducial.root --cats all --doZeroes --ext _Run2 --translateCats cats.json

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
for cmd in "${commands[@]}"; do
    echo "Executing command: $cmd"
    # 执行命令
    eval "$cmd"
    # 检查命令执行状态
    if [ $? -eq 0 ]; then
        echo "Command executed successfully"
    else
        echo "Error executing command: $cmd"
        # 可选择在错误发生时终止执行
        # exit 1
    fi
done

# # 在此处添加代码以等待 Condor 作业完成
# echo "Waiting for Condor job to finish..."
# while true; do
#     if condor_q | grep -q "0 jobs"; then
#         echo "All Condor jobs finished"
#         break
#     else
#         sleep 60  # 每隔60秒轮询一次
#     fi
# done

# commands_1=(
#     "python CollectFits.py --inputJson inputs.json --mode mu_fiducial"
#     "plot1DScan.py runFits_mu_fiducial/profile1D_syst_r_higgs_in.root --y-cut 20 --y-max 20 --output r_fiducial_fixed_statsyst --POI r_higgs_in --translate ../Plots/pois_mu.json --main-label 'Expected' --main-color 1 --others runFits_mu_fiducial/profile1D_statonly_r_higgs_in.root:'Stat only':2 --logo-sub 'Work in Progress'"
# )


# # 遍历命令列表，依次执行每个命令
# for cmd in "${commands_1[@]}"; do
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

echo "All commands executed"