This file describes the content of this folder. The folder includes all R-scripts to apply the workflow developed in the Master's thesis of Nathalie Heiniger with the title:
Identifying Anthropogenic Feeding Sites from GPS Tracking Data: A Case Study for Red Kites (Milvus milvus) in Western Switzerland. (October 2020) 

The workflow includes 5 R-scripts ordered according to their order of application.

Workflow

RMD files
1. 01_msc_nheiniger_torun_datacleaning:
	File includes:
	- download of data from movebank.org
	- data cleaning process (red kite GPS tracking data and additional data for forest and building layers)	
	- result of this R-script is loaded into 02_msc_nheiniger_torun_analysis_KDE.Rmd
	
2. 02_msc_nheiniger_torun_analysis_KDE:
	File includes:
	- preparation of information about private feeders
	- calculation of 50% home range polygons with 'adehabitat'
	- intersection of the polygons with a building layer (prepared in 01_msc_nheiniger_torun_datacleaning)
	- voluntary evaluation of the kernel density estimation approach
	- result of this R-script is loaded into 03_msc_nheiniger_torun_analysis_revisitation.Rmd

3. 03_msc_nheiniger_torun_analysis_revisitation.Rmd:
	File includes:
	- preparation of data for 'recurse' package
	- calculation of revisitation measures
	- result of this R-script is loaded into 05_msc_nheiniger_torun_analysis_GLMM.Rmd

(4. 04_msc_nheiniger_torun_analysis_modelselection.Rmd
	File includes:
	- model selection process
	- code does not need to be run! The 'best' model evaluated in the master's thesis can be loaded from the file "model.rdata")

5. 05_msc_nheiniger_torun_analysis_GLMM.Rmd:
	File includes:
	- preparation of data 
	- predictions based on the 'best' model (file: model.rdata)
	- analysis of model performance
	- evaluation of classification results (true positives, false positives, false negatives)



Folders

- data
This folder contains additional data needed for the error-free implementation of the code files above.
Further information can be found in the README.txt file in the data folder.


- RCode_originals
This folder contains the original, not revised or restructured R-code developed and used for the Master's thesis and 
the file "rk_all_studyarea_NOTforest_daytimeonly_relocs_2018.Rdata" which contains the clean data set used for the Master's thesis analyses.
