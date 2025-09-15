# -*- coding: utf-8 -*-
import ROOT
import os
import glob

# 定义年份列表
years = ["2016preVFP", "2016postVFP", "2017", "2018", "2022preEE", "2022postEE", "2023preBPix", "2023postBPix"]

cat = "ggH0"

# 初始化总统计变量
# 原始统计
total_ele_entries = 0
total_ele_sum = 0.0
total_mu_entries = 0
total_mu_sum = 0.0

# 质量窗口统计 [120, 130]
total_ele_in_window = 0
total_ele_sum_in_window = 0.0
total_mu_in_window = 0
total_mu_sum_in_window = 0.0

# 遍历所有年份
for year in years:
    # 构建文件夹路径
    dir_path = "input_new/signal/signal_{}/".format(year)
    
    # 查找所有包含"M125"的ROOT文件
    root_files = glob.glob(os.path.join(dir_path, "*M125*.root"))
    
    for file_path in root_files:
        print "\nProcessing file:", file_path
        
        # 打开ROOT文件
        tfile = ROOT.TFile.Open(file_path)
        if not tfile or tfile.IsZombie():
            print "  ERROR: Failed to open file"
            continue
        
        # 获取workspace
        ws = tfile.Get("tagsDumper/cms_hgg_13TeV")
        if not ws:
            print "  ERROR: Workspace not found"
            tfile.Close()
            continue
        
        # 定义要读取的数据集名称
        datasets = {
            "ggH_125_ele_13TeV_{}".format(cat): "ele",
            "ggH_125_mu_13TeV_{}".format(cat): "mu",
            "VBF_125_ele_13TeV_{}".format(cat): "ele",
            "VBF_125_mu_13TeV_{}".format(cat): "mu",
        }
        
        # 读取并统计数据集
        for ds_name, ds_type in datasets.items():
            dataset = ws.data(ds_name)
            if not dataset:
                print "  WARNING: Dataset '{}' not found".format(ds_name)
                continue
            
            # 获取原始统计信息
            n_entries = dataset.numEntries()
            sum_entries = dataset.sumEntries()
            
            # print "  Dataset: {} (Type: {})".format(ds_name, ds_type)
            # print "    Original numEntries =", n_entries
            # print "    Original sumEntries =", sum_entries
            
            # 获取质量变量
            mass_var = dataset.get().find("CMS_hgg_mass")
            if not mass_var:
                print "    WARNING: CMS_hgg_mass variable not found in dataset"
                continue
                
            # 创建质量窗口选择器
            mass_var.setRange("mass_window", 120, 130)
            cut_expr = ROOT.RooFit.CutRange("mass_window")
            
            # 创建满足质量窗口的子数据集
            reduced_ds = dataset.reduce(cut_expr)
            if not reduced_ds:
                print "    ERROR: Failed to create reduced dataset for mass window"
                continue
                
            # 获取质量窗口内的统计信息
            n_in_window = reduced_ds.numEntries()
            sum_in_window = reduced_ds.sumEntries()
            
            # print "    In mass window [120, 130]:"
            # print "      numEntries =", n_in_window
            # print "      sumEntries =", sum_in_window
            
            # 累加到总统计
            if ds_type == "ele":
                total_ele_entries += n_entries
                total_ele_sum += sum_entries
                total_ele_in_window += n_in_window
                total_ele_sum_in_window += sum_in_window
            elif ds_type == "mu":
                total_mu_entries += n_entries
                total_mu_sum += sum_entries
                total_mu_in_window += n_in_window
                total_mu_sum_in_window += sum_in_window
            
            # 清理临时对象
            del reduced_ds
        
        # 关闭文件
        tfile.Close()

# 打印总统计结果
print "\n" + "="*70
print "TOTAL STATISTICS ACROSS ALL YEARS:"
print "  Electron datasets:"
print "    Original:"
print "      Total numEntries =", total_ele_entries
print "      Total sumEntries =", total_ele_sum
print "    In mass window [120, 130]:"
print "      Total numEntries =", total_ele_in_window
print "      Total sumEntries =", total_ele_sum_in_window
print "    Fraction in window: {:.2%}".format(float(total_ele_in_window)/total_ele_entries if total_ele_entries else 0)
print "    Weighted fraction: {:.2%}".format(total_ele_sum_in_window/total_ele_sum if total_ele_sum else 0)
print ""
print "  Muon datasets:"
print "    Original:"
print "      Total numEntries =", total_mu_entries
print "      Total sumEntries =", total_mu_sum
print "    In mass window [120, 130]:"
print "      Total numEntries =", total_mu_in_window
print "      Total sumEntries =", total_mu_sum_in_window
print "    Fraction in window: {:.2%}".format(float(total_mu_in_window)/total_mu_entries if total_mu_entries else 0)
print "    Weighted fraction: {:.2%}".format(total_mu_sum_in_window/total_mu_sum if total_mu_sum else 0)
print "="*70
print ""
print "  Total datasets:"
print "    Original:"
print "      Total numEntries =", total_ele_entries + total_mu_entries
print "      Total sumEntries =", total_ele_sum + total_mu_sum
print "    In mass window [120, 130]:"
print "      Total numEntries =", total_ele_in_window + total_mu_in_window
print "      Total sumEntries =", total_ele_sum_in_window + total_mu_sum_in_window
print "    Fraction in window: {:.2%}".format(float(total_ele_in_window + total_mu_in_window)/float(total_ele_entries + total_mu_entries) if (total_ele_entries + total_mu_entries) else 0)
print "    Weighted fraction: {:.2%}".format(total_ele_sum_in_window + total_mu_sum_in_window/float(total_ele_sum + total_mu_sum) if (total_ele_sum + total_mu_sum) else 0)
print "="*70
print "\nProcessing completed."