function OutFile = ReCluster_Cluster(TClust2, RunName)
%% Ben Ireland, June 2026, Cluster multi-variate datasets
% School of Earth Sciences, University of Bristol
%
% Script to recluster clusters from previous clusterings:
%
% Methods in these scripts are described in:
% "Clustering global volcano deformation
% datasets: insights and limitations for analogue signals"
% Submitted to JVGR
%
% Copyright: Ben Ireland, 2026
% V1.0    
    %% Colormaps
    load colorblind_colormap.mat
    load UniqueColors.mat
    
    IBMColor = [100 143 255];
    IBMColor = IBMColor./255;
    
    % Paul Tol's 'muted' colourblind colourmap
    %('www.nceas.ucsb.edu/sites/default/files/2022-06/Colorblind%20Safe%20Color%20Schemes.pdf')
    Colorblind2 = [...
        221 221 221;...
        46 37 133;...
        51 117 56;...
        93 203 236;...
        148 203 236;...
        220 205 125;...
        194 106 119;...
        159 74 150;...
        126 041 084];
    Colorblind2 = Colorblind2./255;
    
    %% Recluster cluster
    RelativeStats=0;
    Save=1;
    
    DistMethod = 'cosine';
    LinkMethod = 'average';
    NumReClust2 = 5;
    
    %% Prep data
    labelsClust2 = TClust2.Volcano;
    TClust2Vals = TClust2(:,2:6);
    % transform and normalise
    TClust2Vals = table2array(TClust2Vals);
    TClust2Vals = asinh(TClust2Vals);
    for i = 1:size(TClust2Vals,2)
        MaxVal2(i)=max(abs(TClust2Vals(:,i)));
        TClust2Vals(:,i) = TClust2Vals(:,i)./MaxVal2(i);
    end
    %% Do clustering
    distances=pdist(TClust2Vals,DistMethod);
    links=linkage(distances, LinkMethod);
    if Save==1
        Folder = pwd;
        mkdir('Re_Clustered_Results')
        Z = links;
        X = TClust2Vals;
        OutFile = [Folder,'/Re_Clustered_Results/ZX_Data_',num2str(RunName),'.mat'];
        save(OutFile,'Z','X');
    end
    Reclust2=cluster(links,'maxclust',NumReClust2);
    %cophenetic correlation coefficient (the closer to 1 the better)
    cClust2=cophenet(links, distances);
    %inconsistency coefficient (higher = more distinct clusters)
    iClust2=inconsistent(links,10);
    maxiClust2=max(iClust2(:,4));
    meaniClust2=mean(iClust2(:,4));
    variClust2=var(iClust2(:,4));
    stdiClust2=std(iClust2(:,4));
        
    cutoff_for_colours = median([links(end-NumReClust2+1,3) links(end-NumReClust2+2, 3)]);
    
    %% Plot Dendrogram
    fig = figure('WindowStyle','docked');
    t1 = tiledlayout(1,1,'TileSpacing','none','Padding','tight');
    nexttile
    [H2, ~, P2, colors,theGroups] = dendrogram2(links, 0, 'Labels',[], 'orientation','left','ColorThreshold',cutoff_for_colours);
    hold on
    Cmap = colorTriplets(2:end,:);
    NumMultClusts2 = recolorDendrogramClasses(H2, Reclust2, labelsClust2, Cmap, 1);
    
    % Add dummy legend
    for j = 1:NumMultClusts2
        dendLgd2(j) = plot(nan,nan,'-','Color',Cmap(j,:),'LineWidth',2,'DisplayName',['Cluster ', num2str(j)]);
    end
    if NumMultClusts2~=NumReClust2
        dendLgd2(end+1) = plot(nan,nan,'-','Color','k','LineWidth',2,'DisplayName','Outlier');
    end
    dendLgd2(end+1) = xline(cutoff_for_colours,'--','LabelHorizontalAlignment','left','LabelVerticalAlignment','top','LabelOrientation','horizontal','DisplayName','Cut-off','LineWidth',2);
    lgdDend = legend(dendLgd2,'Location','southwest'); 
    
    set(gca,'fontsize',10);
    set(H2,'LineWidth',2);
    set(gca,'XTickLabelRotation',0);
    set(gca,'YTickLabelRotation',0);
    xlabel('Linkage height')
    ylabel('Volcano')
    %title('Re-clustering cluster 2 - Dendrogram')
    
    
    %% Re-order clusters
    [GC, GR] = groupcounts(Reclust2);
    OutlierCount = find(GC==1);
    OutlierNum = GR(OutlierCount);
    
    % Create sorted labels
    XLab = string(1:NumReClust2);
    
    %OutlierIdx = find(Reclust2==OutlierNum);
    for j = 1:length(OutlierNum)
        if OutlierNum(j)==1
            XLab = [strcat('O_',num2str(j)), XLab];
        else
            XLab = [XLab(1:OutlierNum(j)-1), strcat('O_',num2str(j)), XLab(OutlierNum(j):end)];
        end
    end
    XLab = XLab(1:NumReClust2);
    [XLabNew, sortedIdx] = sort(XLab,'asc');
    XLabOld = XLab;
    XLab = XLabNew;
    
    % Apply this mapping to clusters
    clustI = Reclust2;
    % Change index of groups and data to match alphabetical XTickLabels
    clustI2 = zeros(size(clustI));
    for j = 1:length(XLab)
        idxThisC = clustI == sortedIdx(j);
        ClustCounts(j) = sum(idxThisC);
        clustI2(idxThisC) = j;
        OldClust = j;
        NewClust = sortedIdx(j);
    end
    clustI = clustI2;
    
    % Work out data rows of outlier clusters
    IdxOutlier = ClustCounts == 1;
    SingleClusts = find(IdxOutlier);
    SingleClustOrder = sort(SingleClusts,'ascend'); % To make appending below simpler
    
    col = colorTriplets(2:end,:); % To avoid black as the first colour
    if ~isempty(SingleClusts)
        for j = 1:length(SingleClusts)
            IdxOutlierData(j) = find(clustI == SingleClusts(j));
    
            % Add black to color triplets for outliers
            col = [col(1:SingleClustOrder(j)-1,:); [0,0,0]; col(SingleClustOrder(j):end,:)];
        end
    end
    % Cut colour array to correct size for number of clusters
    col = col(1:NumReClust2,:);
    
    
    
    %% Plot silhouette scores
    fig = figure('WindowStyle','docked');   
    %title('Re-clustering cluster 2 - Silhouette scores')
    hold on
    [SVals, SFig] = silhouette(TClust2Vals,clustI,DistMethod)
    view([90 -90])
    set(gca,'YDir','normal')
    box on
    
    % Get axes
    ax = SFig.Children;
    %XLab = string(1:NumReClust2);
    ax.YTickLabel = XLab;
    ax.XAxisLocation = 'top';
    SilBar = ax.Children;
    XD = SilBar.XData;
    YD = SilBar.YData;
    
    % Colour bars
    nBars = numel(YD);
    C = repmat(SilBar.FaceColor, nBars, 1);
    
    % Find difference between xdata to find place for XLines
    Gaps = diff(YD);
    GapLocs = find(isnan(YD));
    OutlierNum = 0;
    for k = 1:((length(GapLocs)./2)-1)
        idx = 2*k;
        idx2 = 2*(k+1);
        YLineLoc = mean(GapLocs(idx-1:idx));
        yline(YLineLoc,'--');
    
        Sz = (GapLocs(idx2) - GapLocs(idx)) + 1;
        C(GapLocs(idx):GapLocs(idx2),:) = repmat(col(k,:),Sz,1);
    end
    SilBar.CData = C;
    SilBar.FaceColor = 'flat';
    
    %% Plot cluster stats
    % Plot cluster stats vs my classification and GVP categories
    TStats = TClust2(:,8:end);
    figure('WindowStyle','Docked')
    TC2 = tiledlayout(ceil(width(TStats)./3),3,"TileSpacing","compact","Padding","compact");
    %title(TC2,'Re-clustering cluster 2 - Cluster Stats')
    for i = 1:width(TStats)
        AllClasses = table2array(unique(TClust2(:,i+7)));
        % Get events of each class in each cluster
        for k = 1:NumReClust2
            Classes2 = TClust2{:,i+7};
            ClustIdx = find(Reclust2==k);
            MyClass{k} = Classes2(ClustIdx);
            SizeClust(k) = length(ClustIdx);
            for j = 1:length(AllClasses)
                try
                    NumClass(k,j) = nnz(strcmp(MyClass{k},AllClasses{j}));
                catch
                    NumClass(k,j) = nnz(strcmp(string(MyClass{k}),string(AllClasses(j))));
                end
            end
            if RelativeStats ==1
                NumClass(k,:) = (NumClass(k,:)./SizeClust(k)).*100;
            end
        end
    
        Colorblind2 = [Colorblind2; 0 0 0];
        % Stacked bar chart
        nexttile
        bx = bar(NumClass,'stacked','FaceColor','flat');
        % Change colours to be colorblind friendly
        for x = 1:length(bx)
            bx(x).CData = Colorblind2(x,:);
        end
        if RelativeStats ==1
            ylabel('Percentage of events')
            hold on
            for k = 1:NumReClust2
                text(k,5,0,['n = ',num2str(SizeClust(k))],Color='k',HorizontalAlignment='center',FontWeight='bold',BackgroundColor='w');
            end
        else
            ylabel('Number of events')
        end
        xticklabels(XLab)
        xlabel('Cluster number')
        try
            lgd_bx = legend(bx,AllClasses,'Location','northeast','AutoUpdate','Off');
        catch
            lgd_bx = legend(bx,string(AllClasses),'Location','northeast','AutoUpdate','Off');
        end
        clear NumClass
    end
end
