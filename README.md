The Ntuplizer folder has the details of the electron skimmer used.

The Analyzer folder has the code to create and fill the histograms and store them in a root file.

#Codes to modify

 HGCNtupleVariables.h : Update branches according to the input file
 
 AnalyzeHGCMuons.cc : Contains main analysis code in the main function
 
 AnalyzeHGCMuons.h : Header file for above where histograms are booked and finally saved in the destructor ( ~AnalyzeHGCMuons() )

#For using the analyzer
 
 To run: ./analyzeHGCMuons runlist.txt Plots.root 1100 1100

#Parameters

 runlist.txt : .txt file with the list of files to run over

 Plots.root : Name of the output root file where the created histograms are saved

 The last two arguments are just placeholders.

 The Plotting subfolder in Analyzer has the plotting scripts used. 
 
 Note : The plotting scripts just read the root file with the histograms and then redesign them according to need. The fitting is also done using the plotting scripts.
 
 For using the plotting scripts: Change the name of the rootfiles in the files argument.
 Update the input parameters like histogram name, rebinning , xlabel, ylabel etc. in the respective .C file.
 
 Run the plotting script on lxplus as follows : 
 
 root -b -q generate1Dplot.C 
