# Config file: options for background fitting

backgroundScriptCfg = {
  
  # Setup
  # 'inputWS':'/eos/user/m/mingtao/workspace/zgamma/CMSSW_10_2_13/src/FlashggFinalFit_HZGamma/Inputdata/outputs_rzou_run3_finalfit/Background_workspace/ws/output_Data_all.root', # location of 'output_Data_Run2.root' file
  # 'inputWS':'/eos/user/m/mingtao/workspace/zgamma/CMSSW_10_2_13/src/FlashggFinalFit_HZGamma/Inputdata/input_jiehan/background/ws/output_Data_all.root', # location of 'output_Data_Run2.root' file
  # 'inputWS':'/eos/user/m/mingtao/workspace/zgamma/CMSSW_10_2_13/src/FlashggFinalFit_HZGamma/Inputdata/input_25Sep/background/ws/output_Data_all.root', # location of 'output_Data_Run2.root' file
  'inputWS':'/eos/user/m/mingtao/workspace/zgamma/CMSSW_10_2_13/src/FlashggFinalFit_HZGamma/Inputdata/input_25Dec/background/ws/output_Data_all.root', # location of 'output_Data_Run2.root' file
  # 'inputWS':'/eos/user/m/mingtao/workspace/zgamma/CMSSW_10_2_13/src/FlashggFinalFit_HZGamma/sync_comb/output_Data_all.root', # location of 'output_Data_Run2.root' file
  ## 'ggH0,ggH1,ggH2,ggH3,VBF0,VBF1,VBF2,VBF3' 'Incl0,Incl1,Incl2,Incl3,Incl4,Incl5,Incl6,Incl7'
  # 'cats':'Incl0,Incl1,Incl2,Incl3,Incl4,Incl5,Incl6,Incl7', # auto: automatically inferred from input ws
  'cats':'ggH1', # auto: automatically inferred from input ws
  'catOffset':1, # add offset to category numbers (useful for categories from different allData.root files)  
  ## 'jiehan' 'rzou_run3' 'final'
  'ext':'26Jan_sync', # extension to add to output directory (e.g. final, new)
  'year':'combined', # Use combined when merging all years in category (for plots)

  # Job submission options
  'batch':'local', # [condor,SGE,IC,local]
  'queue':'hep.q' # for condor e.g. microcentury
  
}