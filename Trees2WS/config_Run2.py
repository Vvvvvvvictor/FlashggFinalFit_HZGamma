# Input config file for running trees2ws

trees2wsCfg = {
  # Name of RooDirectory storing input tree
  'inputTreeDir':'DiphotonTree',

  # Variables to be added to dataframe: use wildcard * for common strings
  'mainVars':["CMS_hgg_mass","weight","weight_central","dZ","*Up","*Down", "*central"],
  'dataVars':["CMS_hgg_mass","weight"], # Vars to be added for data
  'stxsVar':'',
  'notagVars':[], # Vars to add to NOTAG RooDataset
  'systematicsVars':["CMS_hgg_mass","weight"], # Variables to add to sytematic RooDataHists
  # 'theoryWeightContainers':{'weight_LHEPdf': 101, 'weight_LHEScale': 9},
  'theoryWeightContainers':{},

  # List of systematics: use string YEAR for year-dependent systematics
  'systematics':["Smearing", "Scale", "FNUF", "Material", "JER", "JES", "METJes", "METUnclustered", "MuonPt"],
  #'systematics': [''],

  # Analysis categories: python list of cats or use 'auto' to extract from input tree
  'cats':'auto'
}
