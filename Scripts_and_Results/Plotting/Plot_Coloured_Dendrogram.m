function Plot_Coloured_Dendrogram(Data)
%% Ben Ireland, June 2026, Cluster multi-variate datasets
% School of Earth Sciences, University of Bristol
%
% Function to plot a dendrogram coloured by feature values
%
% Methods in these scripts are described in:
% "Clustering global volcano deformation
% datasets: insights and limitations for analogue signals"
% Submitted to JVGR
%
% Copyright: Ben Ireland, 2026
% V1.0

load('UniqueColors.mat');
N = 9; % Number of clusters
C = colorTriplets; % Colour triplets
L = {'Duration','Velocity','Area','Aspect Ratio','Mean Depth'}; % Variable labels

for k = 1:length(Data)
    load(strcat(Data(k).folder,'/',Data(k).name));
    Combo = extractBetween(Data(k).name,'Data_','.mat');
    Combo = Combo{:};
    colored_dendrogram_features2(Z, X, L, N, C, Combo, 0);
end