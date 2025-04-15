#!/usr/bin/env python
import os

def rooArgSetToList(argset): ## taken from Andrea Marini's great repo here: https://github.com/amarini/rfwsutils/blob/master/wsutils.py#L300-L313
    """creates a python list with the contents of argset (which should be a RooArgSet)"""
    it = argset.createIterator()

    retval = []
    while True:
        obj = it.Next()

        if obj == None:
            break

        retval.append(obj)

    return retval

def raiseMultiError(lax=False):
    raise RuntimeError('Found more than one multipdf here - please create a workspace with just one for these bias studies. You can use "combineCards.py Datacard.txt --ic cat_name" for this)')

def raiseFailError(itoy, lax=False):
    text = 'some fits have failed, wrong quantile for toy number %g'%itoy
    if not lax: raise RuntimeError('ERROR %s'%text)
    else: print 'WARNING %s'%text

def shortName(name):
    return name.split('_')[-1]

def toyName(name, split=None, jobName=""):
    retval = '/eos/home-j/jiehan/finalfit_102X/CMSSW_10_2_13/src/flashggFinalFit/Combine/Checks/BiasStudy/%s/BiasToys/biasStudy_%s_toys.root' % (jobName, name)
    if split is not None: 
        split = int(split)
        retval = retval.replace(name,'%s_split%g'%(name,split))
    return retval

def fitName(name, split=None, jobName=""):
    retval = '/eos/home-j/jiehan/finalfit_102X/CMSSW_10_2_13/src/flashggFinalFit/Combine/Checks/BiasStudy/%s/BiasFits/biasStudy_%s_fits.root' % (jobName, name)
    if split is not None: 
        split = int(split)
        retval = retval.replace(name,'%s_split%g'%(name,split))
    return retval

def plotName(name, jobName=""):
    return '/eos/home-j/jiehan/finalfit_102X/CMSSW_10_2_13/src/flashggFinalFit/Combine/Checks/BiasStudy/%s/BiasPlots/biasStudy_%s_pulls' % (jobName, name)

def run(cmd, dry=False):
   print cmd
   if not dry: os.system(cmd)

def writeCondorScript(cmds, jobName, dry=False, type='toys'):
    condorFile = open('%s_%s.sub' % (jobName, type), 'w')
    condorFile.write('universe = vanilla\n')
    condorFile.write('Executable = %s_%s.sh\n'%(jobName, type))
    condorFile.write('Arguments = $(ProcId)\n')
    condorFile.write('Output = %s.$(ProcId).out\n'%jobName)
    condorFile.write('Error = %s.$(ProcId).err\n'%jobName)
    condorFile.write('Log = %s.$(ProcId).log\n'%jobName)
    condorFile.write('on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)\n')
    condorFile.write('periodic_release =  (NumJobStarts < 3) && ((CurrentTime - EnteredCurrentStatus) > 600)\n')
    condorFile.write('+JobFlavor = "workday"\n')
    condorFile.write('Queue %s\n' % len(cmds))
    condorFile.close()

    scriptFile = open('%s_%s.sh' % (jobName, type), 'w')
    scriptFile.write('#!/bin/bash\n')
    scriptFile.write('ulimit -s unlimited\n')
    scriptFile.write('set -e\n')
    scriptFile.write('cd %s\n' % os.getcwd())
    scriptFile.write('export SCRAM_ARCH=slc7_amd64_gcc700\n')
    scriptFile.write('source /cvmfs/cms.cern.ch/cmsset_default.sh\n')
    scriptFile.write('eval `scramv1 runtime -sh`\n\n')
    scriptFile.write('cd BiasStudy/%s/%s\n' %(jobName, 'BiasToys' if type == 'toys' else 'BiasFits'))

    for i, cmd in enumerate(cmds):
        scriptFile.write('if [ $1 -eq %s ]; then\n' % i)
        for c in cmd.split('echo ;'):
            scriptFile.write('  %s\n' % c)
        scriptFile.write('fi\n')
    scriptFile.close()

    os.system('chmod +x %s_%s.sh' % (jobName, type))
    if not dry:
        os.system('condor_submit %s_%s.sub' % (jobName, type))