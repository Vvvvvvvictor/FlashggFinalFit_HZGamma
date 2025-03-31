# Config file: options for signal fitting

_year = '2018'

signalScriptCfg = {
  
  # Setup
  'inputWSDir':'/eos/home-j/jiehan/root/input_finalfit/signal_%s'%_year,
  'procs':'auto', # if auto: inferred automatically from filenames
  'cats':'auto', # if auto: inferred automatically from (0) workspace
  'ext':'fiducial_%s'%_year,
  'analysis':'fiducialAnalysis', # To specify which replacement dataset mapping (defined in ./python/replacementMap.py)
  'year':'%s'%_year, # Use 'combined' if merging all years: not recommended
  'massPoints':'120,125,130',
  'flavours': 'ele,mu',

  #Photon shape systematics  
  'scales':'Scale,MuonPt', # separate nuisance per year
  'scalesCorr':'Material,FNUF', # correlated across years
  'scalesGlobal':'', #'NonLinearity,Geant4', # affect all processes equally, correlated across years
  'smears':'Smearing', # separate nuisance per year

  # Job submission options
  'batch':'local', # ['condor','SGE','IC','local']
  'queue':'hep.q'
  #'batch':'condor', # ['condor','SGE','IC','local']
  #'queue':'espresso',

}
