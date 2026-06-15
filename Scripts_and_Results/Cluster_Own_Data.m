%% Ben Ireland, June 2026, Cluster multi-variate datasets
% School of Earth Sciences, University of Bristol
%
% Script to cluster your own dataset (volcanological or not):
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

%% Clustering dataset
% Example dataset with 100 events with random values across 6
% characteristics
Catalogue = 'Example_Input_Dataset.xlsx';
%Catalogue = 'Paper_Dataset_Example.xlsx';

T = readtable(Catalogue,Sheet='Data'); % Data
TInfo = readtable(Catalogue,Sheet='Information'); % Background info to categorise events
T_GroundTruth = readtable(Catalogue,Sheet='Ground_Truth'); % Ground truth

% If you do not have background info or ground truth, uncomment the below lines:
%TInfo = [];
%T_GroundTruth = [];



%% Clustering options
% Transform values before clustering?
% 0 = Don't transform
% 1 = log of absolute values
% 2 = log values and preserve sign
% 3 = Inverse hyperbolic sine (IHS) transform
Lg = 3;

% Normalisation method
% 0 = Don't normalise
% 1 = Z-Score
% 2 = Log normalisation
% 3 = Min/Max normalisation
Normalisation = 3;

% Remove outliers in dataset? (Outside of 2.5% and 97.% percentiles)
NoOutliers=0;

% Give columns to use for clustering and comparison
ClustCols = 2:width(T); % All columns (assuming first column is name)
T = T(:,ClustCols);
if ~isempty(TInfo)
    CompCols = 2:width(TInfo); % Change if you only want to compare to some columns (not all)
    TInfo = TInfo(:,CompCols);
end
if ~isempty(T_GroundTruth)
    GroundCols = 2:width(T_GroundTruth); % Change if you only want to compare to some columns (not all)
    T_GroundTruth = T_GroundTruth(:,GroundCols);
end


%% Validate clustering approach
% Validation-specific options
% Number of clusters to test for validation step
RunName = 'Test_My_Own_Data';
%RunName = 'Test_Paper_Example';
MinClust = 2; % Min num of clusters
MaxClust = 25; % Max num of clusters
ClustGap = 1; % Step in cluster size range
ClustNums = MinClust:ClustGap:MaxClust;
OutFolder = [pwd,'/Example_Results/Approach_Validation']; 
% Folder to save results in - files in there can be fed into the 
% Plot_Coloured_Dendrogram function to visualise

% Tests:
% Distance metrics: Euclidean, Cityblock, Cosine
% Linkage criteria: Single, Average, Complete, Ward's
Validate_Clustering_Approach(T,Lg,Normalisation,NoOutliers,ClustNums,OutFolder,RunName);

%% Do clustering with chosen clustering approach
% Options
numclust = 9; % Can specify multiple if you want to test multiple numbers of clusters at once e.g. [9, 10]
DistMethod = 'cosine';
LinkMethod = 'average';
SaveClusters = 1; % Save clusters to put into coloured dendrogram
SaveClustStats = 1; % Save cluster stats to .csv
OutFolder2 = [pwd,'/Example_Results/Clustering']; 
labels = string(T{:,1}); % Labels/names for each event - for the dendrogram

% Apply chosen clustering method to your dataset, with additional figures
% for analysis
OutFile = Test_Chosen_Clustering_Approach(T_GroundTruth,T,TInfo,RunName,SaveClusters,SaveClustStats,OutFolder2,NoOutliers,Lg,Normalisation,LinkMethod,DistMethod,numclust,labels);

%% Plot coloured dendrogram of clustered results
Data = dir(OutFile);
L = T.Properties.VariableNames;
C = load('UniqueColors.mat').colorTriplets;
for k = 1:length(numclust)
    N = numclust(k);
    Plot_Coloured_DendrogramV2(Data,L,C,N);
end
