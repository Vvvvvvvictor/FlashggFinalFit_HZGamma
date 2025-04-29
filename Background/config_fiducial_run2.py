# Config file: options for signal fitting

backgroundScriptCfg = {
  
  # Setup
  'inputWS':'/eos/home-j/jiehan/root/input_finalfit/background/ws/output_Data_all.root', # location of 'output_Data_Run2.root' file
  'cats':'VBF2', # auto: automatically inferred from input ws
  'catOffset':2, # add offset to category numbers (useful for categories from different allData.root files)  
  'ext':'fiducialAnalysis', # extension to add to output directory
  'year':'combined', # Use combined when merging all years in category (for plots)

  # Job submission options
  'batch':'local', # [condor,SGE,IC,local]
  'queue':'hep.q' # for condor e.g. microcentury
  
}