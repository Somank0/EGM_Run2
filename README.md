> [!NOTE]
> All the instructions here are to perform the DRN energy regression on the AOD data.

## To create the EDProducer : gsfElectronDRN :

1. Go to CMSSW_X_Y_Z/src

2. git cms-addpkg PhysicsTools/PatAlgos

3. copy the ***patElectronDRNCorrector_cfi.py*** and place it in the *PhysicaTools/PatAlgis/python/slimming* 

> [!TIP]
> For the studies, **CMSSW_13_3_3** was used, but any CMSSW_version >= 13_3_X should work.

## Inference script

The executable *DRN_job_withFileCheck.sh* is used to run the inference. It calls the file DRN_reg_final_cfg.py 

To run inference files of a dataset :

```
./DRN_job_withFileCheck.sh [folder_name] [dataset_name] [start_file_index] [end_file_index]
```

To run inference on first 'n' files of a dataset :

```
./DRN_job_withFileCheck.sh [folder_name] [dataset_name] 0 n
```

This was done to batch the jobs for condor.

A sample condor script is written in the file : ***EGM_condorjob_pfThreshTL_pedTL.sub***


The Ntuplizer folder has the details of the electron skimmer used.

> [!IMPORTANT]
> The DRN_job_withFileCheck.sh file calls the nTuplizer. Hence the nTuplizer must be created with the name ***Electron_RefinedRecHit_NTuplizer***. 

The Analyzer folder has the code to create and fill the histograms and store them in a root file.

## Codes to modify

 HGCNtupleVariables.h : Update branches according to the input file
 
 AnalyzeHGCMuons.cc : Contains main analysis code in the main function
 
 AnalyzeHGCMuons.h : Header file for above where histograms are booked and finally saved in the destructor ( ~AnalyzeHGCMuons() )

## For using the analyzer
 
 To run: ./analyzeHGCMuons runlist.txt Plots.root 1100 1100

### Parameters

 runlist.txt : .txt file with the list of files to run over

 Plots.root : Name of the output root file where the created histograms are saved

 The last two arguments are just placeholders.

 The Plotting subfolder in Analyzer has the plotting scripts used. 
 
 Note : The plotting scripts just read the root file with the histograms and then redesign them according to need. The fitting is also done using the plotting scripts.
 
 For using the plotting scripts: Change the name of the rootfiles in the files argument.
 Update the input parameters like histogram name, rebinning , xlabel, ylabel etc. in the respective .C file.
 
 Run the plotting script on lxplus as follows : 
 
 root -b -q generate1Dplot.C 


## To create a nTuplizer :

1. Go to the CMSSW_X_Y_Z/src

2. Run the following commands :

```
cmsenv

mkedanlzr Electron_RefinedRecHit_NTuplizer

cd Electron_RefinedRecHit_NTuplizer/plugins
```
3. Replace the .h, .cc and .xml files in the plugins folder with the .h, .cc and .xml file in the NTuplizer folder of this git respository.

4. Run the following command :
```
scram b -j32
```

This creates the Ntuplizer plugin to be called by DRN_reg_final_cfg.py.

