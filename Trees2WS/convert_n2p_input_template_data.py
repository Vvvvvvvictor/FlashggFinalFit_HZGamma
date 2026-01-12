# -*- coding: UTF-8 -*-
import ROOT
import glob
import os
from array import array

def process_data_background():
    # 定义BDT分界线
    ggF_bdt_bins = [0.57, 0.83, 0.94]
    VBF_bdt_bins = [0.48, 0.81, 0.91]
    
    # 创建输出目录
    output_dir = "Inputdata/input_25Dec/templates"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # 创建输出文件
    output_file = ROOT.TFile(os.path.join(output_dir, "sstest_templates_data.root"), "RECREATE")
    
    # 初始化直方图字典
    histograms = {}
    histograms_sb = {}
    
    # 定义直方图类别
    categories = ["ggH0", "ggH1", "ggH2", "ggH3", "VBF0", "VBF1", "VBF2", "VBF3"]
    
    # 为每个类别创建直方图
    for cat in categories:
        hist_name = "data_full_%s" % cat
        hist_name_sb = "data_%s" % cat
        histograms[cat] = ROOT.TH1F(hist_name, hist_name, 340, 95, 180)  # 340个bin，范围95-180 GeV
        histograms_sb[cat] = ROOT.TH1F(hist_name_sb, hist_name_sb, 340, 95, 180)
    
    # 定义要处理的文件模式
    file_patterns = [
        "data_*_output.root"
    ]
    
    # 要排除的文件模式
    exclude_patterns = [
        # 排除的文件
    ]
    
    # 处理ggF和VBF两个目录
    inputfilePath = "/eos/project/h/htozg-dy-privatemc/rzou/bdt/BDT_output_redwood"
    directories = [
        "%s/Output_ggF_rui_redwood_v1_ext_val"%(inputfilePath), 
        "%s/Output_VBF_rui_redwood_v1_ext_val"%(inputfilePath)
        ]
    
    for dir_path in directories:
        print("处理目录: %s" % dir_path)
        
        for file_pattern in file_patterns:
            # 构建完整的文件路径模式
            full_pattern = os.path.join(dir_path, file_pattern)
            files = glob.glob(full_pattern)
            
            # 排除不需要的文件
            filtered_files = []
            for file_path in files:
                exclude = False
                for exclude_pattern in exclude_patterns:
                    exclude_full_pattern = os.path.join(dir_path, exclude_pattern)
                    exclude_files = glob.glob(exclude_full_pattern)
                    if file_path in exclude_files:
                        exclude = True
                        break
                if not exclude:
                    filtered_files.append(file_path)
            
            print("在目录 %s 中找到 %d 个 %s 文件" % (dir_path, len(filtered_files), file_pattern))
            
            for file_path in filtered_files:
                print("处理文件: %s" % file_path)
                try:
                    input_file = ROOT.TFile(file_path, "READ")
                    tree = input_file.Get("outtree")
                    
                    n_entries = tree.GetEntries()
                    
                    # 遍历所有条目
                    for i in range(n_entries):
                        tree.GetEntry(i)
                        
                        # 直接访问branch的值
                        llphoton_m = tree.llphoton_refit_m
                        weight_corr = tree.weight_corr
                        bdt_val = tree.BDT_score
                        
                        # 根据目录确定过程类型
                        if "ggF" in dir_path:
                            process_type = "ggF"
                            bdt_bins = ggF_bdt_bins
                        else:  # VBF目录
                            process_type = "VBF" 
                            bdt_bins = VBF_bdt_bins
                        
                        # 根据BDT分数分类确定类别
                        if process_type == "ggF":
                            # ggF分类
                            if bdt_val < bdt_bins[0]:
                                category = "ggH3"
                            elif bdt_val < bdt_bins[1]:
                                category = "ggH2"
                            elif bdt_val < bdt_bins[2]:
                                category = "ggH1"
                            else:
                                category = "ggH0"
                        else:
                            # VBF分类（根据您的要求：>0.91→VBF3, 0.81-0.91→VBF2, 0.48-0.81→VBF1, <0.48→VBF0）
                            if bdt_val > bdt_bins[2]:  # > 0.91
                                category = "VBF0"
                            elif bdt_val > bdt_bins[1]:  # > 0.81
                                category = "VBF1"
                            elif bdt_val > bdt_bins[0]:  # > 0.48
                                category = "VBF2"
                            else:  # < 0.48
                                category = "VBF3"
                        
                        # 填充直方图
                        histograms[category].Fill(llphoton_m, weight_corr)
                        if ((llphoton_m < 120) or (llphoton_m > 130)):
                            histograms_sb[category].Fill(llphoton_m, weight_corr)
                    
                    input_file.Close()
                            
                except Exception as e:
                    print("处理文件 %s 时出错: %s" % (file_path, e))
    
    # 写入直方图
    print("开始写入直方图...")
    output_file.cd()
    
    for cat in categories:
        hist = histograms[cat]
        if hist.GetEntries() > 0:
            hist.Write()
            print("写入直方图 %s, 包含 %.1f 个加权事件" % (hist.GetName(), hist.Integral()))
        else:
            print("直方图 %s 没有数据，跳过" % cat)

    for cat in categories:
        hist_sb = histograms_sb[cat]
        if hist_sb.GetEntries() > 0:
            hist_sb.Write()
            print("写入直方图 %s, 包含 %.1f 个加权事件" % (hist_sb.GetName(), hist_sb.Integral()))
        else:
            print("直方图 %s 没有数据，跳过" % cat)
   
    output_file.Close()
    print("处理完成！输出文件: %s" % os.path.join(output_dir, "sstest_templates_data.root"))

if __name__ == "__main__":
    process_data_background()