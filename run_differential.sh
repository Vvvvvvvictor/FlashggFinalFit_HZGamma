#!/bin/bash

# 定义要执行的命令列表
commands=(
    # 环境
    "source /cvmfs/cms.cern.ch/cmsset_default.sh"
    "cmsenv"
    "source setup.sh"
    # Tree2WS
    # "cd Trees2WS/"
    # "bash run_differential.sh"
    # "cd /eos/user/c/chpan/input_finalfit/"
    # "cp hh.sh signal_differential_2022postEE"
    # "cd signal_differential_2022postEE"
    # "bash hh.sh"
    # "mv hh.sh ../signal_differential_2022preEE"
    # "cd ../signal_differential_2022preEE"
    # "bash hh.sh"
    # "rm -rf hh.sh"
    # "cd /eos/user/c/chpan/finalfit/CMSSW_10_2_13/src/flashggFinalFit/"
    
    # Signal
    "cd Signal"
    # fTest
    # "python RunSignalScripts.py --inputConfig config_2022preEE_differential.py --mode fTest --modeOpts "--doPlots""
    # "python RunSignalScripts.py --inputConfig config_2022postEE_differential.py --mode fTest --modeOpts "--doPlots""
    # # syst
    # "python RunSignalScripts.py --inputConfig config_2022preEE_differential.py --mode calcPhotonSyst"
    # "python RunSignalScripts.py --inputConfig config_2022postEE_differential.py --mode calcPhotonSyst"
    # signalfit
    "python RunSignalScripts.py --inputConfig config_2022preEE_differential.py --mode signalFit --groupSignalFitJobsByCat --modeOpts "--doPlots""
    "python RunSignalScripts.py --inputConfig config_2022postEE_differential.py --mode signalFit --groupSignalFitJobsByCat --modeOpts "--doPlots""
    # packaged
    "python RunPackager.py --cats best_resolution_0jets,best_resolution_1jets,best_resolution_2jets,best_resolution_3jets,medium_resolution_0jets,medium_resolution_1jets,medium_resolution_2jets,medium_resolution_3jets,worst_resolution_0jets,worst_resolution_1jets,worst_resolution_2jets,worst_resolution_3jets --exts differential_2022preEE,differential_2022postEE --mergeYears --batch local"
    # Background
    # "cd ../Background"
    # "python RunBackgroundScripts.py --inputConfig config_2022_differential.py --mode fTestParallel"
    # Datacard
    # "cd ../Datacard"
    # "python RunYields.py --inputWSDirMap 2022preEE=/eos/user/c/chpan/input_finalfit/signal_differential_2022preEE/,2022postEE=/eos/user/c/chpan/input_finalfit/signal_differential_2022postEE/ --cats auto --procs auto --doSystematics --batch local --ext differential --mergeYears --skipZeroes"
    # "python makeDatacard.py --years 2022preEE,2022postEE --ext differential --prune --pruneThreshold 0.001 --doSystematics --doMCStatUncertainty"
    # # Combine
    # "cd ../Combine"
    # "cp ../Datacard/Datacard.txt ."
    # "mkdir Models"
    # "mkdir Models/signal"
    # "mkdir Models/background"
    # "cp -r ../Signal/outdir_packaged/CMS-HGG_sigfit_packaged_* Models/signal/"
    # "cp -r ../Background/outdir_differentialAnalysis/CMS-HGG_multipdf_* Models/background/"
    # # 改改名字
    # "sed -i -E '/13TeV_bkgshape/s/(13TeV_bkgshape)/2022_\1_norm/g; /pdfindex/s/(resolution)/\1_2022/g' Datacard.txt"
    # # 送到afs去
    # "cp Datacard.txt ~/finalfit/CMSSW_10_2_13/src/flashggFinalFit/Combine/"
    # "cp -r Models ~/finalfit/CMSSW_10_2_13/src/flashggFinalFit/Combine/"
    # "cd ~/finalfit/CMSSW_10_2_13/src/flashggFinalFit/"
    # # 环境
    # "source /cvmfs/cms.cern.ch/cmsset_default.sh"
    # "cmsenv"
    # "source setup.sh"
    # # Combine
    # "cd Combine"
    # "python RunText2Workspace.py --mode mu_differential --batch local"
    # "python RunFits.py --inputJson inputs.json --mode mu_differential --batch condor"

    # 在这里添加更多命令
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

# 在此处添加代码以等待 Condor 作业完成
echo "Waiting for Condor job to finish..."
while true; do
    if condor_q | grep -q "0 jobs"; then
        echo "All Condor jobs finished"
        break
    else
        sleep 60  # 每隔60秒轮询一次
    fi
done

# commands_1=(
#     "python CollectFits.py --inputJson inputs.json --mode mu_differential"
#     "plot1DScan.py runFits_mu_differential/profile1D_syst_r_higgs_in_0jets.root --y-cut 20 --y-max 20 --output plot_differential_0617/r_higgs_in_0jets --POI r_higgs_in_0jets --translate ../Plots/pois_mu.json --main-label "Expected" --main-color 1 --others runFits_mu_differential/profile1D_statonly_r_higgs_in_0jets.root:"Stat only":2 --logo-sub "Work in Progress""
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

# echo "All commands executed"
