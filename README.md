# Final Fits (dev_higgsdnafinalfit)

**NOTE** This is the devlopment branch for using final fits with the output of HiggsDNA(not updated HiggsDNA version).

Welcome to the new Final Fits package. Here lies a a series of scripts which are used to run the final stages of the CMS Hgg analysis: signal modelling, background modelling, datacard creation, final statistical interpretation and final result plots.

Slides from the flashgg tutorial series can be found [here](https://indico.cern.ch/event/963619/contributions/4112177/attachments/2151275/3627204/finalfits_tutorial_201126.pdf)

## Download and setup instructions

```
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_2_13
cd CMSSW_10_2_13/src

# Install the GBRLikelihood package which contains the RooDoubleCBFast implementation
git clone https://github.com/jonathon-langford/HiggsAnalysis.git

# Install Combine as per the documentation here: cms-analysis.github.io/HiggsAnalysis-CombinedLimit/
git clone -b v8.2.0 https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit

# Install Combine Harvester for parallelizing fits
git clone -b 102x https://github.com/cms-analysis/CombineHarvester.git CombineHarvester

# Compile external libraries
cmsenv
scram b -j 9

# Install Flashgg Final Fit packages
git clone -b dev_higgsdnafinalfit https://github.com/cms-analysis/flashggFinalFit.git
cd flashggFinalFit/
```

In every new shell run the following to add `tools/commonTools` and `tools/commonObjects` to your `${PYTHONPATH}`:
```
cmsenv
source setup.sh
```

## Contents
The Finals Fits package contains several subfolders which are used for the following steps:

* Create the Signal Model (see `Signal` dir)
* Create the Background Model (see `Background` dir)
* Generate a Datacard (see `Datacard` dir)
* Running fits with combine (see `Combine` dir)
* Scripts to produce plots (see `Plots` dir)

The signal modelling, background modelling and datacard creation can be ran in parallel. Of course the final fits (`Combine`) requires the output of these three steps. In addition, the scripts in the `Trees2WS` dir are a series of lightweight scripts for converting standard ROOT trees into a RooWorkspace that can be read by the Final Fits package.

Finally, the objects and tools which are common to all subfolders are defined in the `tools` directory. If your input workspaces differ from the flashgg output workspace structure, then you may need to change the options here.

Each of the relevant folders are documented with specific `README.md` files. Some (temporary) instructions can be found in this [google docs](https://docs.google.com/document/d/1NwUrPvOZ2bByaHNqt_Fr6oYcP7icpbw1mPlw_3lHhEE/edit)


# Important illustration of the final fits workflow

## Generate the workspace for fitting
example for the data workspace
```
cd Trees2WS
python trees2ws/trees2ws.py --inputCondfig <data_config_file> --inputTreeFile <specific_format_file> --outputWSDir <output_dir>
python trees2ws/trees2ws.py --inputCondfig <signal_config_file> --inputTreeFile <specific_format_file> --inputMass <mass_point> --productionMode <production_mode> --year <year> --outputWSDir <output_dir> (--doSystematics)
```

## Create the signal model

### Make signal model
```
cd Signal
python RunSignalScripts.py --inputConfig <config_file> --mode signalFit --groupSignalFitJobsByCat --modeOpts "--doPlots --beamspotWidthData 3.5 --beamspotWidthMC 3.7 --useDCB --skipSystematics"
```
*Only use the DCB function to fit*

### Package the signal model
```
python RunPackager.py --cats <cats-to-combine> --exts <ext-document-in-config-file> --mergeYears --batch local --outputExt <output_ext>
```

### Create the signal model workspace
```
python RunPlotter.py --procs all --cats all --years <years-to-package> --ext <output_ext> 
```

### Create the signal pdf in cats for spurious signal test
```
python scripts/combineSignalPdf.py --cat <cat>
```

## Create the background model
```
python RunBackgroundScripts.py --inputConfig <config_file> --mode fTestParallel --jobOpts “--blindFit”
```
*blind fit for unblind period*

**Check the paths in fTest.cpp/SpurialSignalTest before run fTest**
* path for spurious signal test MC shape: `TFile* fbkg = TFile::Open(Form("/eos/user/j/jiehan/root/input_finalfit/templates/template_%s.root", runPeriod.Data()));`, if you don't have one, please generate it with `/HiggsZGammaAna/synchronization_script/create_direct_templates.py` or `HiggsZGammaAna/SSTest/Generate_template.py`
* path for signal combined pdf path: `TString signalFileName = Form("/eos/user/j/jiehan/finalfit_102X/CMSSW_10_2_13/src/flashggFinalFit/Signal/outdir_combinedPDFs/CMS-HGG_combinedPDFs_%s.root", cat.Data());`

## Make the datacard
```
cd Datacard
python RunYields.py --inputWSDirMap <signal_mc_ws_file> --cats <category_to_process/auto> --ext <yield_ext> ----mergeYears --skipZeroes (--doSystematics)
python mackDatacard.py --years <{year}={path-to-signal-model-ws},...> --prune --pruneThreshold 0.001 --mergeYears --output <output_name> (--doSystematics)
```
* with `--prune` option, the datacard will be pruned to remove the categories with less than 0.1% production mode

## Run Bias study
```
cd Combine
python RunText2Workspace.py --ext <datacard_ext>  --mode <mode-in-model.py> --batch local
python RunBiasStudy.py -j <job-name>  -d <path-to-datacard> -t -n <ntoys> (--dryRun)

cd Check
python RunBiasStudy.py -j <job-name>  -d <path-to-datacard> -f -n <ntoys>  -c "--cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints   --rMin -100 --rMax 100  --freezeParameters MH" (--dryRun)
python RunBiasStudy.py -j <job-name>  -d <path-to-datacard> -p --gaussianFit
```

**-t: generate toys; -f: fit the toys with envelope; -p: plot the test result**

## Run the final fit
```
python RunFits.py --inputJson <input_json_file> --ext <datacard_ext> --mode <mode-in-model.py> --batch condor/local
```

**If do fitting in each mu point**
```
python CollectFits.py --inputJson <input_json_file> --ext <datacard_ext> --mode <mode-in-model.py>
plot1DScan.py <final_fit_output> --y-cut <line-threshold> --y_max <figure-threshold> --output <syst-output-name> --POI r --translate <poi-json> --main-label 'Expected' --main-color 1 --others <statonly-output-name>:'Stat only':2 --logo-sub 'Work in Progress'
```

## Make the S+B plot
```
python makeSplusBModelPlot.py --inputWSFile <datacard-workspace> --cats all --doZeroes --ext <output_ext> --translateCats <cat-json>