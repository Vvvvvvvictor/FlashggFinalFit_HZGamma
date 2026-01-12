# -*- coding: utf-8 -*-
import ROOT
import os
import sys

def deep_clone_tree(tree, directory=None):
    """安全地深克隆一个树"""
    if not tree:
        return None
    
    # 先克隆树结构
    cloned_tree = tree.CloneTree(0)  # 0表示只克隆结构
    
    # 设置目录
    if directory:
        cloned_tree.SetDirectory(directory)
    
    # 确保树被正确设置
    cloned_tree.SetName(tree.GetName())
    
    return cloned_tree

def merge_data_files():
    # 年份列表
    years = ["2016preVFP", "2016postVFP", "2017", "2018", 
             "2022preEE", "2022postEE", "2023preBPix", "2023postBPix"]
    
    # 定义要处理的树名称
    categories = ["VBF", "ggH"]
    subcategories = ["0", "1", "2", "3"]
    flavors = ["ele", "mu"]
    
    # 输出文件
    output_dir = "/eos/user/m/mingtao/workspace/zgamma/CMSSW_10_2_13/src/FlashggFinalFit_HZGamma/Inputdata/input_25Dec/fitting_bkg/Data"
    output_filename = "output_Data_all.root"
    
    # 创建输出目录
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    output_path = os.path.join(output_dir, output_filename)
    
    # 如果文件已存在，删除它
    if os.path.exists(output_path):
        print("Removing existing output file...")
        os.remove(output_path)
    
    # 打开输出文件
    output_file = ROOT.TFile(output_path, "CREATE")
    
    # 创建DiphotonTree目录
    output_dir_diphoton = output_file.mkdir("DiphotonTree")
    
    # 用于存储合并树的对象
    merged_trees = {}
    
    # 为每个可能的树名创建占位符
    all_tree_names = []
    for cat in categories:
        for subcat in subcategories:
            for flav in flavors:
                tree_name = "Data_13TeV_{cat}{subcat}_{flav}".format(
                    cat=cat, subcat=subcat, flav=flav
                )
                all_tree_names.append(tree_name)
                merged_trees[tree_name] = None
    
    print("Processing files...")
    
    # 处理每个年份的文件
    for year_idx, year in enumerate(years):
        input_filename = "/eos/user/m/mingtao/workspace/zgamma/CMSSW_10_2_13/src/FlashggFinalFit_HZGamma/Inputdata/input_25Dec/fitting_bkg/Data_{year}/output_Data_{year}.root".format(year=year)
        
        print("\n[{idx}/{total}] Processing {filename}".format(
            idx=year_idx+1, total=len(years), filename=input_filename))
        
        if not os.path.exists(input_filename):
            print("  Warning: File does not exist. Skipping.")
            continue
        
        # 打开输入文件
        input_file = ROOT.TFile(input_filename, "READ")
        
        if input_file.IsZombie():
            print("  Error: Cannot open file. Skipping.")
            input_file.Close()
            continue
        
        # 获取DiphotonTree目录
        diphoton_dir = input_file.Get("DiphotonTree")
        if not diphoton_dir:
            print("  Warning: DiphotonTree directory not found. Skipping.")
            input_file.Close()
            continue
        
        # 遍历所有可能的树
        for tree_name in all_tree_names:
            # 获取输入树
            input_tree = diphoton_dir.Get(tree_name)
            if not input_tree:
                continue
            
            n_entries = input_tree.GetEntries()
            if n_entries == 0:
                continue
            
            print("  Found {tree_name}: {entries} entries".format(
                tree_name=tree_name, entries=n_entries))
            
            # 检查是否需要创建输出树
            if merged_trees[tree_name] is None:
                print("    Creating output tree for {tree_name}".format(tree_name=tree_name))
                # 切换到输出目录
                output_dir_diphoton.cd()
                # 深克隆树结构
                merged_trees[tree_name] = deep_clone_tree(input_tree, output_dir_diphoton)
                if not merged_trees[tree_name]:
                    print("    ERROR: Failed to create output tree")
                    continue
            
            # 获取输出树
            output_tree = merged_trees[tree_name]
            
            # 手动复制所有条目
            entries_added = 0
            try:
                print("    Copying entries...")
                
                # 使用手动逐个复制（最可靠的方法）
                total_entries = input_tree.GetEntries()
                
                # 分批处理大文件
                batch_size = 10000
                for start_idx in range(0, total_entries, batch_size):
                    end_idx = min(start_idx + batch_size, total_entries)
                    
                    for i in range(start_idx, end_idx):
                        input_tree.GetEntry(i)
                        output_tree.Fill()
                        entries_added += 1
                    
                    # 显示进度
                    if end_idx % 100000 == 0 or end_idx == total_entries:
                        print("      Processed {current}/{total} entries ({percent:.1f}%)".format(
                            current=end_idx, total=total_entries, 
                            percent=100.0 * end_idx / total_entries))
                
                print("    Added {added} entries from {year}".format(
                    added=entries_added, year=year))
                
            except Exception as e:
                print("    ERROR during copy: {error}".format(error=str(e)))
                import traceback
                traceback.print_exc()
        
        # 关闭输入文件
        input_file.Close()
        
        # 显示当前统计
        print("  Current statistics after {year}:".format(year=year))
        for tree_name in all_tree_names:
            tree = merged_trees[tree_name]
            if tree:
                n = tree.GetEntries()
                if n > 0:
                    print("    {tree_name}: {entries} entries".format(
                        tree_name=tree_name, entries=n))
    
    # 写入所有合并的树
    print("\n\nWriting merged trees to output file...")
    output_dir_diphoton.cd()
    
    total_entries = 0
    trees_written = 0
    
    for tree_name in all_tree_names:
        tree = merged_trees[tree_name]
        if tree is not None:
            n_entries = tree.GetEntries()
            if n_entries > 0:
                # 确保在正确的目录中
                output_dir_diphoton.cd()
                
                # 写入树
                tree.Write(tree_name, ROOT.TObject.kOverwrite)
                total_entries += n_entries
                trees_written += 1
                print("  Written {tree_name}: {entries} entries".format(
                    tree_name=tree_name, entries=n_entries))
            else:
                print("  Skipping {tree_name}: 0 entries".format(tree_name=tree_name))
        else:
            print("  No data for {tree_name}".format(tree_name=tree_name))
    
    # 保存并关闭输出文件
    output_file.Write()
    output_file.Close()
    
    print("\nMerge completed!")
    print("Output file: {path}".format(path=output_path))
    print("Total entries: {entries}".format(entries=total_entries))
    print("Trees written: {count}".format(count=trees_written))
    
    # 验证输出文件
    print("\nVerifying output file...")
    verify_file = ROOT.TFile(output_path, "READ")
    if not verify_file.IsZombie():
        diphoton_dir_out = verify_file.Get("DiphotonTree")
        if diphoton_dir_out:
            print("Trees in output file:")
            total_verified = 0
            trees_verified = 0
            
            for key in diphoton_dir_out.GetListOfKeys():
                obj = key.ReadObj()
                if isinstance(obj, ROOT.TTree):
                    n = obj.GetEntries()
                    print("  {name}: {entries} entries".format(
                        name=obj.GetName(), entries=n))
                    total_verified += n
                    trees_verified += 1
            
            print("\nVerification summary:")
            print("  Trees found: {trees}".format(trees=trees_verified))
            print("  Total entries: {entries}".format(entries=total_verified))
            
            if total_entries != total_verified:
                print("WARNING: Entry count mismatch!")
                print("  Expected: {expected}".format(expected=total_entries))
                print("  Actual: {actual}".format(actual=total_verified))
                print("  Difference: {diff}".format(diff=total_entries - total_verified))
        else:
            print("Error: DiphotonTree directory not found in output")
    else:
        print("Error: Cannot open output file for verification")
    
    verify_file.Close()

if __name__ == "__main__":
    # 启用批处理模式
    ROOT.gROOT.SetBatch(True)
    
    # 增加缓冲区大小以提高性能
    ROOT.TTree.SetMaxTreeSize(500 * 1024 * 1024 * 1024)  # 500GB
    
    # 禁用一些可能导致问题的特性
    ROOT.gEnv.SetValue("TFile.AsyncPrefetching", 0)
    ROOT.gEnv.SetValue("TFile.AsyncReading", 0)
    
    # 增加堆栈大小
    import resource
    resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))
    
    # 设置递归限制
    import sys
    sys.setrecursionlimit(10000)
    
    merge_data_files()