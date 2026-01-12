# -*- coding: UTF-8 -*-
import ROOT
import glob
import os
from array import array

def process_all_files():
    # 定义BDT分界线
    ggF_bdt_bins = [0.57, 0.83, 0.94]
    VBF_bdt_bins = [0.48, 0.81, 0.91]
    
    # 创建输出目录
    output_dir = "Inputdata/input_25Sep/fitting_signal/all_M125_all"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # 创建输出文件
    output_file = ROOT.TFile(os.path.join(output_dir, "output_all_M125.root"), "RECREATE")
    
    # 创建TDirectoryFile
    diphoton_dir = output_file.mkdir("DiphotonTree")
    diphoton_dir.cd()
    
    # 初始化数据存储字典
    data_dict = {}
    
    # 定义TTree名称模板
    ggf_ele_trees = ["all_125_ele_13TeV_ggH%d" % i for i in range(4)]
    ggf_mu_trees = ["all_125_mu_13TeV_ggH%d" % i for i in range(4)]
    vbf_ele_trees = ["all_125_ele_13TeV_VBF%d" % i for i in range(4)]
    vbf_mu_trees = ["all_125_mu_13TeV_VBF%d" % i for i in range(4)]
    
    all_tree_names = ggf_ele_trees + ggf_mu_trees + vbf_ele_trees + vbf_mu_trees
    
    # 为每个TTree初始化数据列表
    for name in all_tree_names:
        data_dict[name] = {
            "CMS_hgg_mass": [],
            "weight": [],
            "bdt_score": []
        }
    
    # 定义要处理的目录和文件模式
    inputfilePath = "/eos/project/h/htozg-dy-privatemc/rzou/bdt/BDT_output_redwood"
    directories = [
        ("%s/Output_ggF_rui_redwood_v1_ext_val"%(inputfilePath), "GGF_*_output.root", "ggF"),
        ("%s/Output_ggF_rui_redwood_v1_ext_val"%(inputfilePath), "VBF_*_output.root", "ggF"), 
        ("%s/Output_VBF_rui_redwood_v1_ext_val"%(inputfilePath), "GGF_*_output.root", "VBF"),
        ("%s/Output_VBF_rui_redwood_v1_ext_val"%(inputfilePath), "VBF_*_output.root", "VBF")
    ]
    
    for dir_path, file_pattern, process_type in directories:
        # 构建完整的文件路径模式
        full_pattern = os.path.join(dir_path, file_pattern)
        files = glob.glob(full_pattern)
        
        print("在目录 %s 中找到 %d 个 %s 文件" % (dir_path, len(files), file_pattern))
        
        for file_path in files:
            print("处理文件: %s (类型: %s)" % (file_path, process_type))
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
                    nel = getattr(tree, 'nel', 0)
                    nmu = getattr(tree, 'nmu', 0)
                    
                    # 根据nel和nmu确定通道
                    if nel == 0 and nmu == 2:
                        channel = "mu"
                    elif nel == 2 and nmu == 0:
                        channel = "ele"
                    else:
                        print("Warning: 异常的事例 nel=%d, nmu=%d，跳过此事例" % (nel, nmu))
                        continue
                    
                    # 根据处理类型选择BDT分界线
                    if process_type == "ggF":
                        bdt_bins = ggF_bdt_bins
                        process_prefix = "ggH"
                    else:  # VBF
                        bdt_bins = VBF_bdt_bins
                        process_prefix = "VBF"
                    
                    # 根据BDT分数分类确定bin
                    if process_type == "ggF":
                        # ggF: 正常顺序
                        if bdt_val < bdt_bins[0]:
                            bdt_bin = 0
                        elif bdt_val < bdt_bins[1]:
                            bdt_bin = 1
                        elif bdt_val < bdt_bins[2]:
                            bdt_bin = 2
                        else:
                            bdt_bin = 3
                    else:
                        # VBF: 反向顺序（根据您之前的描述）
                        if bdt_val > bdt_bins[2]:  # > 0.91
                            bdt_bin = 3
                        elif bdt_val > bdt_bins[1]:  # > 0.81
                            bdt_bin = 2
                        elif bdt_val > bdt_bins[0]:  # > 0.48
                            bdt_bin = 1
                        else:
                            bdt_bin = 0
                    
                    # 确定TTree名称
                    tree_name = "all_125_%s_13TeV_%s%d" % (channel, process_prefix, bdt_bin)
                    
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
    for tree_name in all_tree_names:
        n_events = len(data_dict[tree_name]["CMS_hgg_mass"])
        print("处理 %s, 包含 %d 个事件" % (tree_name, n_events))
        
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
    print("处理完成！输出文件: %s" % os.path.join(output_dir, "output_all_M125.root"))

if __name__ == "__main__":
    process_all_files()