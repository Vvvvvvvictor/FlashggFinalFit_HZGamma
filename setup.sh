# Add tools dir to PYTHONPATH
PS1=' [\W] $ '
eval `scramv1 runtime -sh`

export PYTHONPATH=$PYTHONPATH:${CMSSW_BASE}/src/FlashggFinalFit_HZGamma/tools
