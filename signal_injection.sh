#!/bin/bash

# 初始化
clear
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsenv
source $(pwd)/setup.sh
mainPath=$(pwd)
cd ${mainPath}/Background
# make clean; make -j 16;

# 设置一些prefix
cat="ggH2"
prefix="final_corr"
output_file="outdir_${prefix}/bkgfTest-Data/signal_inject.txt"
function="PowerLawStepxGau3"
# sed -i "s#final_2#${prefix}#g" config_fiducial_run2.py

# 清空或创建输出文件
echo "$cat $function" > "$output_file"

# 循环处理每个i值
for i in 0 1 2 3 5 10; do

    file="outdir_${prefix}/bkgfTest-Data/all_${cat}_${function}_${i}xsig.txt"
    # rm $file
    # signal injection test
    echo "Start signal injection test"
    # python RunBackgroundScripts.py --inputConfig config_fiducial_run2.py --mode fTestParallel --jobOpts "--nsig_in $i"  #“--blindFit”

    # 检查文件是否存在
    if [ ! -f "$file" ]; then
        echo "$output_file does not exist"
        echo "i = $i: *" >> "$output_file"
        continue
    fi
    
    # 提取目标数值
    value=$(grep -oP "sig \(unweighted\): \K-?\d+(\.\d+)?" "$file" 2>/dev/null | head -n1 | xargs -I {} printf "%.2f" {})
    error=$(grep -oP "sig \(unweighted\): -?\d+(\.\d+)? \+\/- \K-?\d+(\.\d+)?" "$file" 2>/dev/null | head -n1 | xargs -I {} printf "%.2f" {})
    # value=$(grep -m1 "sig (unweighted):" "$file" 2>/dev/null | awk '{printf "%.2f", $3}')
    # error=$(grep -m1 "sig (unweighted):" "$file" 2>/dev/null | awk '{printf "%.2f", $5}')

    # 处理未找到值的情况
    if [ -z "$value" ]; then
        echo "i = $i: *" >> "$output_file"
    else
        echo "i = $i: ${value} +/- ${error}" >> "$output_file"  # 写入i和Ns的值
    fi
done

echo "Results have been saved in $output_file"