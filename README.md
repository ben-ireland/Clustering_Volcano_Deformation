Scripts for clustering multi-variate datasets and analysing results in MATLAB using agglomerative hierarchical clustering.
Author: Ben Ireland, University of Bristol
Copyright: Ben Ireland, 2026
Date: June 2026
Version: 1.0
Contact: Ben Ireland, b.ireland@bristol.ac.uk 

Reference:
Ireland, B., Biggs, J., and Anantrasirichai, N. (2026) Clustering global volcano deformation datasets: insights and limitations for analogue signals. Submitted

For full details of files and scripts in this repository, please see the README.pdf file.

Quick start:

Note: Before running RecreateClusteringResults.m or Cluster_Own_Data.m, ensure the current working directory is “*/Scripts_and_Results” (use pwd command to check, and cd command to change if needed).

1.	I want to reproduce the figures and clustering from the manuscript

  a.	In MATLAB, open the file “Scripts_and_Results  RecreateClusteringResults.m”
  
  b.	Run RecreateClusteringResults.m script. Plots will appear as in the paper, and results can be found in the Cluster_Validation_Results, ClusterStats, Paper_Clustering_Results, and Re_Clustered_Results folders. See Files and file structure section for more detail.


3.	I want to try different agglomerative hierarchical clustering approaches on the same dataset

  a.	In MATLAB, open the file “Scripts_and_Results  Cluster_Own_Data.m”

  
  b.	Ensure the line “Catalogue = 'Paper_Dataset_Example.xlsx';” is uncommented (Line 22).

  
  c.	On lines 86 and 87, change DistMethod and LinkMethod to your desired distance metric and linkage criteria. See https://uk.mathworks.com/help/stats/pdist.html and mathworks.com/help/stats/linkage.html for all options
  d.	Run Cluster_Own_Data.m script. Visualisation validation and clustering plots will appear and results can be found in the Example_Results folder. See Files and file structure section for more detail.

5.	I want to see the input volcano deformation catalogue
  a.	Download the “Deformation_Catalogues  Deformation_Catalogue.xlsx” file. There are alternative versions in the same folder showing only the events used in the clustering, and removing any interpolated values, see Files and file structure section for more detail.
  b.	The volcano deformation events in these catalogues come from common events between the catalogues of Biggs and Pritchard (2017), and Ebmeier et al. (2018).

6.	I want to apply the same clustering approaches to my own multi-variate dataset
  a.	Prepare your dataset in the format shown in the file “Deformation_Catalogues  Example_Dataset  Example_Input_Dataset.xlsx”.
  b.	In MATLAB, open the file “Scripts_and_Results  Cluster_Own_Data.m”
  c.	On Line 21, specify your input dataset .xlsx file prepared in the previous step.
  d.	If you do not have any background information and/or ground truth, uncomment lines 29 and 30.
  e.	Specify your clustering options, run name, and output folder in lines 34-74 following descriptions in the file.
  f.	Run Cluster_Own_Data.m script. Visualisation validation and clustering plots will appear and results can be found in the Example_Results folder. See Files and file structure section for more detail.

7.	I want to do something else
  a.	Please see Files and file structure and individual MATLAB scripts and functions for more detail. If you have a specific request, you can contact the author at the address at the top of the this document.
