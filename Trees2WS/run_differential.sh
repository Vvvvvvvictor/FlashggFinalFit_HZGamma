#!/bin/bash

# 定义参数组合
channels=("GluGluHtoGG" "ttHtoGG" "VBFHtoGG" "VHtoGG")
production_modes=("ggh" "tth" "vbf" "vh")
posts=("postEE" "preEE")
# posts=("postEE")
mass_points=("120" "125" "130")

# 循环遍历所有参数组合
for channel_index in "${!channels[@]}"; do
  channel="${channels[channel_index]}"
  production_mode="${production_modes[channel_index]}"

  for post in "${posts[@]}"; do
    for mass_point in "${mass_points[@]}"; do
      # 构建目录路径
      directory="/eos/user/c/chpan/root_differential_0602/${channel}_M-${mass_point}_${post}/"

      # 检查目录是否存在
      if [ -d "$directory" ]; then
        # 获取目录下所有文件
        files=("$directory"/*)

        # 循环遍历所有文件
        for file in "${files[@]}"; do
          # 构建命令
          command="python trees2ws.py --inputConfig config_2022.py --inputTreeFile $file --inputMass ${mass_point} --productionMode ${production_mode} --year 2022${post} --doSystematics --outputWSDir /eos/user/c/chpan/input_finalfit/signal_differential_2022${post}/"

          # 打印并执行命令
          echo "Executing: $command"
          eval $command

          # 换行
          echo
        done
      else
        echo "Directory not found: $directory"
      fi
    done
  done
done
