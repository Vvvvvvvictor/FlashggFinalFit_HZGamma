# Config file: options for signal fitting

_year = 'all'

signalScriptCfg = {
  
  # Setup
  'inputWSDir':'/eos/home-m/mingtao/workspace/zgamma/CMSSW_10_2_13/src/FlashggFinalFit_HZGamma/Inputdata/input_25Dec/signal/signal_all',
  'procs':'all', # if auto: inferred automatically from filenames
  'cats':'ggH0,ggH1,ggH2,ggH3,VBF0,VBF1,VBF2,VBF3', # if auto: inferred automatically from (0) workspace
  # 'cats':'ggH0,ggH1,ggH2,ggH3',
  'ext':'fiducial_%s'%_year,
  'analysis':'fiducialAnalysis', # To specify which replacement dataset mapping (defined in ./python/replacementMap.py)
  'year':'%s'%_year, # Use 'combined' if merging all years: not recommended
  'massPoints':'120,125,130',
  'flavours': 'ele,mu', # ele,mu,all

  # Shape systematics  
  'scales':'', # separate nuisance per year
  'scalesCorr':'', # correlated across years
  'scalesGlobal':'', #'NonLinearity,Geant4', # affect all processes equally, correlated across years
  'smears':'', # separate nuisance per year

  # Job submission options
  'batch':'local', # ['condor','SGE','IC','local']
  'queue':'hep.q'
  #'batch':'condor', # ['condor','SGE','IC','local']
  #'queue':'espresso',

}
