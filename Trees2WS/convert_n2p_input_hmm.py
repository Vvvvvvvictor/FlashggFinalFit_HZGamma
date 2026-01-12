# -*- coding: UTF-8 -*-
import ROOT
import os

def convert_roodatahist():
    # 创建输出目录
    output_dir = "../Inputdata/input_25Sep/bkg_hmm"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # 创建输出文件
    output_file = ROOT.TFile(os.path.join(output_dir, "output_all_M125.root"), "RECREATE")
    
    # 创建tagsDumper目录
    tags_dumper_dir = output_file.mkdir("tagsDumper")
    tags_dumper_dir.cd()
    
    # 创建主workspace
    main_workspace = ROOT.RooWorkspace("cms_hgg_13TeV", "cms_hgg_13TeV")
    
    # 定义输入输出类别映射
    input_to_output_map = {
        "ggf1": "ggH3",
        "ggf2": "ggH2", 
        "ggf3": "ggH1",
        "ggf4": "ggH0",
        "vbf1": "VBF3",
        "vbf2": "VBF2",
        "vbf3": "VBF1", 
        "vbf4": "VBF0"
    }
    
    # 输入文件
    input_file = ROOT.TFile("/afs/cern.ch/user/m/mioshiro/public/test_datacard11_rawdata.root", "READ")
    
    # 首先在输出workspace中创建统一的变量
    output_var_name = "CMS_hgg_mass"
    output_var = ROOT.RooRealVar(output_var_name, output_var_name, 100, 180)  # 范围可以根据需要调整
    getattr(main_workspace, 'import')(output_var)
    
    # 处理每个类别
    for input_cat, output_cat in input_to_output_map.items():
        print("处理类别: %s -> %s" % (input_cat, output_cat))
        
        # 输入workspace名称
        input_ws_name = "WS_Htomm_cat_%s" % input_cat
        input_datahist_name = "mcdata_Htomm_cat_%s_nominal" % input_cat
        input_var_name = "mllg_cat_%s" % input_cat
        
        # 输出RooDataHist名称
        output_datahist_name = "all_125_all_13TeV_%s" % output_cat
        
        try:
            # 从输入文件获取workspace
            input_ws = input_file.Get(input_ws_name)
            if not input_ws:
                print("错误: 未找到workspace %s" % input_ws_name)
                continue
            
            # 从workspace中获取RooDataHist和变量
            input_datahist = input_ws.data(input_datahist_name)
            input_var = input_ws.var(input_var_name)
            
            if not input_datahist:
                print("错误: 未找到RooDataHist %s" % input_datahist_name)
                continue
            
            if not input_var:
                print("错误: 未找到变量 %s" % input_var_name)
                continue
            
            # 创建新的RooDataHist，使用统一的变量
            # 首先获取原始RooDataHist的bin信息和权重
            hist = input_datahist.createHistogram("temp_hist", input_var)
            n_bins = hist.GetNbinsX()
            
            # 创建新的RooDataHist
            output_datahist = ROOT.RooDataHist(output_datahist_name, output_datahist_name, 
                                             ROOT.RooArgSet(output_var))
            
            # 填充数据到新的RooDataHist
            for i in range(1, n_bins + 1):  # ROOT hist bins start from 1
                bin_center = hist.GetBinCenter(i)
                weight = hist.GetBinContent(i)
                error = hist.GetBinError(i)
                
                if weight > 0:
                    output_var.setVal(bin_center)
                    output_datahist.add(ROOT.RooArgSet(output_var), weight)
                    bin_index = output_datahist.getHist().FindBin(bin_center)
                    output_datahist.getHist().SetBinError(bin_index, error)
            
            # 导入到输出workspace
            getattr(main_workspace, 'import')(output_datahist)
            
            print("成功创建 %s, 包含 %d 个bin, 总权重 %.3f" % 
                  (output_datahist_name, n_bins, output_datahist.sumEntries()))
            
        except Exception as e:
            print("处理类别 %s 时出错: %s" % (input_cat, e))
    
    # 写入输出workspace
    print("写入输出workspace...")
    tags_dumper_dir.cd()
    main_workspace.Write()
    
    # 关闭文件
    input_file.Close()
    output_file.Close()
    
    print("处理完成！输出文件: %s" % os.path.join(output_dir, "output_all_M125.root"))

if __name__ == "__main__":
    convert_roodatahist()