%% Ben Ireland, June 2026, Cluster multi-variate datasets
% School of Earth Sciences, University of Bristol
%
% Script to reproduce figures and results from the manuscript below:
%
% Methods in these scripts are described in:
% "Clustering global volcano deformation
% datasets: insights and limitations for analogue signals"
% Submitted to JVGR
%
% Copyright: Ben Ireland, 2026
% V1.0

% Note: may need to run multiple times if the first time you get an "layout has
% insufficient space" error.

clear all; close all; clc;
addpath(genpath(fileparts(pwd)))

%% Create Upset plot
T = readtable('Deformation_Catalogue.xlsx',Sheet='Full_Catlogue');
T2 = readtable('Deformation_Catalogue_No_Interpolation.xlsx');
T3 = readtable('Deformation_Catalogue.xlsx',Sheet='ColsForClustering');
TInfo = readtable('Deformation_Catalogue.xlsx',Sheet='VolcInfo');
GVP = readtable('GVP_Volcano_List_Holocene.xlsx');
GVP2 = readtable('GVP_Volcano_List_Pleistocene.xlsx');
compareInterp = 1;

Make_UpSetPlot_Fig(T,T2,T3,TInfo,GVP,GVP2,compareInterp);

clear T T2 T3 TInfo GVP GVP2 compareInterp

%% Validate clustering approach
T = readtable('Deformation_Catalogue.xlsx',Sheet='ColsForClustering_NoDuplicates');
TInfo = readtable('Deformation_Catalogue.xlsx',Sheet='VolcInfo_NoDuplicates');

OutFolder = ValidateClusteringApproach_Paper(T,TInfo);

clear T TInfo

%% Do clustering with cosine distance and average linkage metric
T = readtable('Deformation_Catalogue.xlsx',Sheet='ColsForClustering_NoDuplicates');
TInfo = readtable('Deformation_Catalogue.xlsx',Sheet='VolcInfo_NoDuplicates');
MyClasses = readtable('Clustered_Events_Only.xlsx',Sheet='Clusters_MyClassification');

OutFile = Do_Average_Cosine_Clustering(MyClasses,T,TInfo);

clear T TInfo MyClasses

%% Plot coloured dendrogram of clustered results
Data = dir(OutFile);
Plot_Coloured_Dendrogram(Data);

clear Data

%% Re-cluster clusters 1 and 4
% Cluster 4
ClustNum = 4;
TClust = readtable('Clustered_Events_Only.xlsx',Sheet='Cluster4_Events');

OutFile = ReCluster_Cluster(TClust, ClustNum);

Data = dir(OutFile);
Plot_Coloured_Dendrogram(Data);

clear Data ClustNum

% Cluster 1
ClustNum = 1;
TClust = readtable('Clustered_Events_Only.xlsx',Sheet='Cluster1_Events');

OutFile = ReCluster_Cluster(TClust, ClustNum);

Data = dir(OutFile);
Plot_Coloured_Dendrogram(Data);