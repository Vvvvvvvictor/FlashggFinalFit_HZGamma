# Add tools dir to PYTHONPATH
PS1=' [\W] $ '
eval `scramv1 runtime -sh`

export PYTHONPATH=$PYTHONPATH:${CMSSW_BASE}/src/flashggFinalFit/tools
export PYTHONPATH=$PYTHONPATH:${CMSSW_BASE}/src/flashggFinalFit/Trees2WS
