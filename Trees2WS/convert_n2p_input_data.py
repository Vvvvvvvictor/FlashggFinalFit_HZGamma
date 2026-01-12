# -*- coding: UTF-8 -*-
import ROOT
import glob
import os

def process_root_files():
    # 定义BDT分界线
    ggF_bdt_bins = [0.57, 0.83, 0.94]
    VBF_bdt_bins = [0.48, 0.81, 0.91]
    
    # 创建输出文件
    output_file = ROOT.TFile("combined_output.root", "RECREATE")
    
    # 创建TDirectoryFile
    diphoton_dir = output_file.mkdir("DiphotonTree")
    diphoton_dir.cd()
    
    # 初始化数据存储字典
    data_dict = {}
    tree_names = [
        "Data_13TeV_ggH0", "Data_13TeV_ggH1", "Data_13TeV_ggH2", "Data_13TeV_ggH3",
        "Data_13TeV_VBF0", "Data_13TeV_VBF1", "Data_13TeV_VBF2", "Data_13TeV_VBF3"
    ]
    
    # 为每个TTree初始化数据列表
    for name in tree_names:
        data_dict[name] = {
            "CMS_hgg_mass": [],
            "weight": [],
            "bdt_score": []
        }

    inputfilePath = "/eos/project/h/htozg-dy-privatemc/rzou/bdt/BDT_output_redwood"

    # 处理ggF文件
    ggF_files = glob.glob("%s/Output_ggF_rui_redwood_v1_ext_val/data_*_output.root"%(inputfilePath))
    print("找到 %d 个ggF文件" % len(ggF_files))
    
    for file_path in ggF_files:
        print("处理ggF文件: %s" % file_path)
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
                
                # 根据BDT分数分类
                if bdt_val < ggF_bdt_bins[0]:
                    tree_name = "Data_13TeV_ggH0"
                elif bdt_val < ggF_bdt_bins[1]:
                    tree_name = "Data_13TeV_ggH1"
                elif bdt_val < ggF_bdt_bins[2]:
                    tree_name = "Data_13TeV_ggH2"
                else:
                    tree_name = "Data_13TeV_ggH3"
                
                # 存储数据
                data_dict[tree_name]["CMS_hgg_mass"].append(llphoton_m)
                data_dict[tree_name]["weight"].append(weight_corr)
                data_dict[tree_name]["bdt_score"].append(bdt_val)
            
            input_file.Close()
                    
        except Exception as e:
            print("处理文件 %s 时出错: %s" % (file_path, e))
    
    # 处理VBF文件
    VBF_files = glob.glob("%s/Output_VBF_rui_redwood_v1_ext_val/data_*_output.root"%(inputfilePath))
    print("找到 %d 个VBF文件" % len(VBF_files))
    
    for file_path in VBF_files:
        print("处理VBF文件: %s" % file_path)
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
                
                # 根据BDT分数分类（注意VBF的分界线顺序不同）
                if bdt_val < VBF_bdt_bins[0]:  # 0.48
                    tree_name = "Data_13TeV_VBF0"
                elif bdt_val < VBF_bdt_bins[1]:  # 0.81
                    tree_name = "Data_13TeV_VBF1"
                elif bdt_val < VBF_bdt_bins[2]:  # 0.91
                    tree_name = "Data_13TeV_VBF2"
                else:
                    tree_name = "Data_13TeV_VBF3"
                
                # 存储数据
                data_dict[tree_name]["CMS_hgg_mass"].append(llphoton_m)
                data_dict[tree_name]["weight"].append(weight_corr)
                data_dict[tree_name]["bdt_score"].append(bdt_val)
            
            input_file.Close()
                    
        except Exception as e:
            print("处理文件 %s 时出错: %s" % (file_path, e))
    
    # 创建并写入TTrees
    print("开始写入输出文件...")
    diphoton_dir.cd()
    
    # 为每个TTree创建并写入数据
    for tree_name in tree_names:
        n_events = len(data_dict[tree_name]["CMS_hgg_mass"])
        print("处理 %s，包含 %d 个事件" % (tree_name, n_events))
        
        if n_events > 0:
            # 创建TTree和变量
            tree = ROOT.TTree(tree_name, tree_name)
            
            # 创建变量（使用数组方式）
            cms_hgg_mass = array('d', [0.0])
            weight_val = array('d', [0.0])
            bdt_score_val = array('d', [0.0])
            
            # 设置branch
            tree.Branch("CMS_hgg_mass", cms_hgg_mass, "CMS_hgg_mass/D")
            tree.Branch("weight", weight_val, "weight/D")
            tree.Branch("bdt_score", bdt_score_val, "bdt_score/D")
            
            # 填充数据
            for i in range(n_events):
                cms_hgg_mass[0] = data_dict[tree_name]["CMS_hgg_mass"][i]
                weight_val[0] = data_dict[tree_name]["weight"][i]
                bdt_score_val[0] = data_dict[tree_name]["bdt_score"][i]
                tree.Fill()
            
            # 写入TTree
            tree.Write()
            print("写入 %s, 包含 %d 个事件" % (tree_name, tree.GetEntries()))
        else:
            print("%s 没有数据，跳过" % tree_name)
    
    output_file.Close()
    print("处理完成！输出文件: combined_output.root")

# 在文件开头添加array导入
from array import array

if __name__ == "__main__":
    process_root_files()