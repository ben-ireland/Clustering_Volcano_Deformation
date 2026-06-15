function OutFile = Test_Chosen_Clustering_Approach(MyClasses,T,TInfo,RunName,SaveClusters,SaveClustStats,OutFolder,NoOutliers,Lg,Normalisation,LinkMethod,DistMethod,numclust,labels)
%% Ben Ireland, June 2026, Cluster multi-variate datasets
% School of Earth Sciences, University of Bristol
%
% Script to test different clustering approaches.
%
% Methods in these scripts are described in:
% "Clustering global volcano deformation
% datasets: insights and limitations for analogue signals"
% Submitted to JVGR
%
% Copyright: Ben Ireland, 2026
% V1.0

%% Colormaps
    colorblind = load('colorblind_colormap.mat').colorblind;
    colorTriplets = load('UniqueColors.mat').colorTriplets;
    
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
    
    
    %% Input and data pre-processing for Clustering
    %Ben Ireland, cluster deformation data using hierarchicial
    %clustering
    
    if ~isempty(TInfo)
        Validation_Figures = 1;
    else
        Validation_Figures = 0;
    end

    if ~isempty(MyClasses)
        Ground_Truth_Figs = 1;
    else
        Ground_Truth_Figs = 0;
    end

    ParamNames = T.Properties.VariableNames;
    NumVars = width(T);
    
    % Creating labels/subsets of string parameters and removing NaNs
    % Remove NaNs
    T = table2array(T);
    NanMask = ~any(isnan(T), 2);
    
    if NoOutliers==1
        OutMask = isoutlier(T,"percentiles",[2.5 97.5],1);
        OutMaskRow = sum(OutMask,2)>0;
    
        TCrop = T(NanMask & ~OutMaskRow,:);
    else
        TCrop = T(NanMask,:);
    end
    
    % Transforming data
    AddRight = 1;
    TCropOrig = TCrop;
    if Lg ==1 % Optionally take log of absolute values
        TCrop = abs(TCrop);
        TCrop = log10(TCrop);
        Lg_name = 'log_';
    elseif Lg==2 % Take log and preserve sign
        Signs = sign(TCrop);
        TCrop = abs(TCrop);
        TCrop = log10(1+TCrop);
        TCrop = TCrop.*Signs;
        Lg_name = 'log_signed_';
    elseif Lg==3 % Inverse hyperbolic sign
        TCrop = asinh(TCrop);
        Lg_name = 'log_IHS_';
    else
        Lg_name = '';
    end
    
    TCropLog = TCrop;
    
    % Normalising data
    if Normalisation ==1
        % ZScore
        mu = mean(TCrop);
        sigma = std(TCrop);
        TCrop = zscore(TCrop);
    elseif Normalisation ==2
        % Take log
        TCrop = abs(TCrop);
        TCrop = log10(TCrop);
    elseif Normalisation ==3
        % Scale between -1 and 1
        for i = 1:size(TCrop,2)
            MaxVal(i)=max(abs(TCrop(:,i)));
            TCrop(:,i) = TCrop(:,i)./MaxVal(i);
        end
    end
    
    %% Visualise distribution of parameters with histograms
    figure('WindowStyle','docked')
    TL = tiledlayout(3,size(TCrop,2),TileIndexing="columnmajor",Padding="tight",TileSpacing="tight");
    count = 'a';
    for k = 1:size(TCrop,2)
        nexttile
        histogram(TCropOrig(:,k),'NumBins',25,'FaceColor',IBMColor)
        subtitle('Original data')
        %title(ParamNames{k})
        box on
        if k==1
            ylabel('Number of events')
        end
    
        title(ParamNames{k})
    
        % Add subplot letter
        drawnow
        xLimits = xlim;
        yLimits = ylim;
        ylim([0 yLimits(2)*1.2]);
        yLimits = [0 yLimits(2)*1.2]; % Adjust xlims so letter doesn't appear on top
        text(xLimits(1) + 0.02*(xLimits(2)-xLimits(1)), ...
            0.98*yLimits(2), ...
            strcat(count,'.'), 'FontWeight', 'bold', 'FontSize', 12, ...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    
        count = char(count+1);
    
        nexttile
        histogram(TCropLog(:,k),'NumBins',25,'FaceColor',IBMColor)
        subtitle('Transformed')
        box on
    
        if k==1
            ylabel('Number of events')
        end
    
        % Add subplot letter
        drawnow
        xLimits = xlim;
        yLimits = ylim;
        ylim([0 yLimits(2)*1.2]);
        yLimits = [0 yLimits(2)*1.2]; % Adjust xlims so letter doesn't appear on top
        text(xLimits(1) + 0.02*(xLimits(2)-xLimits(1)), ...
            0.98*yLimits(2), ...
            strcat(count,'.'), 'FontWeight', 'bold', 'FontSize', 12, ...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    
        count = char(count+1);
    
        nexttile
        histogram(TCrop(:,k),'NumBins',25,'FaceColor',IBMColor)
        subtitle('Transformed and normalised')
        box on
    
        if k==1
            ylabel('Number of events')
        end
        xlabel('Feature Value')
    
        % Add subplot letter
        drawnow
        xLimits = xlim;
        yLimits = ylim;
        ylim([0 yLimits(2)*1.2]);
        yLimits = [0 yLimits(2)*1.2]; 
        text(xLimits(1) + 0.02*(xLimits(2)-xLimits(1)), ...
            0.98*yLimits(2), ...
            strcat(count,'.'), 'FontWeight', 'bold', 'FontSize', 12, ...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    
        count = char(count+1);
    end
    fontsize(gcf, 14, 'points')


    %% Cluster data Hierarchically
    distances=pdist(TCrop,DistMethod);
    links=linkage(distances, LinkMethod);

    if SaveClusters==1
        Z = links;
        X = TCrop;
        mkdir('Paper_Clustering_Results')
        OutFile = [OutFolder,'/ZX_Data_',RunName,'.mat'];
        save(OutFile,'Z','X');
    end

    %% Plot dendrogram
    for k = 1:length(numclust)
        clust(:,k)=cluster(links,'maxclust',numclust(k));

        %cophenetic correlation coefficient (the closer to 1 the better)
        c(k)=cophenet(links, distances);
        %inconsistency coefficient (higher = more distinct clusters)
        I(k,:,:)=inconsistent(links,10);
        maxi(k)=max(I(k,:,4));
        meani(k)=mean(I(k,:,4));
        vari(k)=var(I(k,:,4));
        stdi(k)=std(I(k,:,4));

        cutoff_for_colours(k) = median([links(end-numclust(k)+1,3) links(end-numclust(k)+2, 3)]);

        %Plot Dendrogram
        fig = figure('WindowStyle','docked');
        t1 = tiledlayout(1,1,'TileSpacing','none','Padding','tight');
        nexttile
        [H, ~, P, colors{k},theGroups(:,k)] = dendrogram2(links, 0, 'Labels',labels, 'orientation','left','ColorThreshold',cutoff_for_colours(k));
        hold on

        T = clust(:,k);                 % your cluster() output
        Cmap = colorTriplets(2:end,:);
        NumMultClusts = recolorDendrogramClassesV2(H, Z, T, Cmap, 1,labels);

        % Add dummy legend
        for j = 1:NumMultClusts
            dendLgd(j) = plot(nan,nan,'-','Color',Cmap(j,:),'LineWidth',2,'DisplayName',['Cluster ', num2str(j)]);
        end

        if NumMultClusts~=numclust(k)
            dendLgd(end+1) = plot(nan,nan,'-','Color','k','LineWidth',2,'DisplayName','Outlier');
        end

        dendLgd(end+1) = xline(cutoff_for_colours(k),'--','LabelHorizontalAlignment','left','LabelVerticalAlignment','top','LabelOrientation','horizontal','DisplayName','Cut-off','LineWidth',2);
        lgdDend = legend(dendLgd,'Location','southwest');

        set(gca,'fontsize',10);
        set(H,'LineWidth',2);
        set(gca,'XTickLabelRotation',0);
        set(gca,'YTickLabelRotation',0);
        xlabel('Linkage height')
        ylabel('Event')
        %title("Clusters: " + numclust(k))

        %box on

        %Get dendrogram colours
        linesColor{k} = cell2mat(get(H,'Color')); % get lines color;
        colorList{k} = unique(linesColor{k}, 'rows');
        clustList{k} = unique(clust(:,k));

        counts_clust{k} = histcounts(clust(:,k),numclust(k));

        % Work out the mapping between theGroups and clust (control the
        % colour order in dendrogram function)

        % Get cluster numbers
        len = length(unique(clust(:,k)));
        Clusts = clust(:,k);
        for j = 1:len
            idx = Clusts==j;
            clustNum(j) = sum(Clusts == j);
            if clustNum(j)==1
                Clusts(idx) = 0; % Put in the same format as theGroups
            end
        end

        Vals1 = unique(Clusts);
        Vals2 = unique(theGroups);
        for j = 1:length(unique(Clusts)) % Compare sums from theGroups and Clusts
            ClustSum(j) = sum(Clusts == Vals1(j));
            GroupSum(j) = sum(theGroups == Vals2(j));
        end
    end

    
    %% Compare clusters to background info and ground truth
    % Prepare collections of groups for each background info category
    if Validation_Figures ==1
        T = TInfo;
        nCols = width(T);     
        allCats = {};
        for c = 1:nCols
            % Extract the column as an array
            colData = T{:, c};
    
            % Ensure it's treated as string
            if iscell(colData)
                colData = string(colData);
            end
    
            % Find unique labels
            Labels{c} = unique(colData);
            allCats = [allCats; Labels{c}];   
        end
    
        T_Variables = T.Properties.VariableNames;
    
        %varNames = T_Variables(1:end-2);
        varNames = T_Variables;
        allCats = unique(allCats); 
        numCats = numel(allCats);
    
        C = zeros(numCats, nCols);   % rows = categories, columns = table columns
        for col = 1:nCols
            colData = T.(varNames{col});
            for k = 1:numCats
                C(k, col) = sum(strcmp(colData, allCats{k}));
            end
        end
    
        %% Plot background stats
        figure('WindowStyle','docked');
        tlStats = tiledlayout(nCols,nCols,"Padding","loose","TileSpacing","tight",'TileIndexing','rowmajor');
    
        nexttile([nCols nCols])
        hold on
        for k = 1:nCols
            SB_Data{k} = C(C(:,k)>0,k);
            h{k} = bar(k,SB_Data{k}, 'stacked');
            for a = 1:length(SB_Data{k})
                h{k}(a).FaceColor = Colorblind2(a,:);
            end
        end
        ylabel('Number of events');
        xlim([0.5 nCols+0.5])
        ylim([0 height(TInfo)]);
        yLims = ylim;
        xticks(1:nCols);
        xticklabels(varNames);
    
        yyaxis right
        ylim(yLims)
        TickVals = [0 height(TInfo)./4 height(TInfo)./2 3*(height(TInfo)./4) height(TInfo)];
        yticks(TickVals)
        yticklabels({'0','25','50','75','100'})
        ylabel('Proportion of signals (%)')
        ax = gca;
        ax.YAxis(1).Color = [0 0 0];
        ax.YAxis(2).Color = [0 0 0];
        box on
    
        %% Plot separate legend for background stats (due to large size)
        figure('WindowStyle','docked');
        tlStats2 = tiledlayout(3,nCols,"Padding","compact","TileSpacing","none",'TileIndexing','rowmajor');
    
        for colnum = 1:nCols
            nexttile([3 1])
            col = colnum;
            % Categories present in this column
            present = C(:,col) > 0;
    
            % Create dummy patches (invisible) for legend
            catsForLegend = allCats(present);
            colsForLegend = Colorblind2(1:length(SB_Data{colnum}), :);
    
            for k = 1:sum(present)
                dummy(k) = patch(nan, nan, colsForLegend(k,:));  % dummy patch
            end
    
            % Plot legend
            catsForLegend = cellfun(@(x) [upper(x(1)) x(2:end)], catsForLegend, 'UniformOutput', false);
            legend(catsForLegend, ...
                'Location', 'north', ...
                'Box', 'off');
    
            % Make invisible
            set(gca,'Visible','off')    
        end
    
        %% Plot stacked bar charts of cluster totals for each background category and ground truth
        for k = 1:length(numclust)
            clustI = clust(:,k); 
    
            fig = figure('WindowStyle','docked');
            numRows = nCols/2;
            if Ground_Truth_Figs ==1
                if mod(nCols,2)==0
                    numRows = numRows + 1;
                end
            end
            TL = tiledlayout(numRows,2,'TileSpacing','compact','Padding','tight');

            if length(numclust)>1
                txt2 = title(TL,"Clusters: " + numclust(k));
                txt2.FontSize = 16;
            end
            letter = 'a'; % For subplot numbers
            for m = 1:nCols
                % Initialise results array
                Results = zeros(numclust(k),length(Labels{m}));
    
               % Extract counts for each variable
                varData = T.(m);
                nvals = numel(Labels{m});
    
                % For each cluster, find counts for each variable
                for l = 1:numclust(k) 
                    Vals = varData(clustI==l); % Extract values for param m for cluster l
                    for n = 1:length(Labels{m})
                        Results(l,n) = sum(matches(Vals,Labels{m}(n))); % Get counts for each label option
                    end
                end
    
                % Make cluster and outlier numbers link with those in the
                % dendrogram and violin plots
                XLab = string(1:numclust);
                OutlierIdx = find(sum(Results')==1);
                for j = 1:length(OutlierIdx)
                   XLab = [XLab(1:OutlierIdx(j)-1), strcat('O_',num2str(j)), XLab(OutlierIdx(j):end)];
                end
                XLab = XLab(1:numclust);
                [XLabNew, sortedIdx] = sort(XLab,'asc');
                XLabOld = XLab;
                XLab = XLabNew;
    
                % plot stacked bar
                nexttile
                b = bar(Results(sortedIdx,:),'stacked','FaceColor','flat');
    
                % Change colours to be colorblind friendly
                for o = 1:length(b)
                    b(o).CData = Colorblind2(o,:);
                end
    
                % Ensure all labels start with captial letter
                Labels2{m} = cellfun(@(x) [upper(x(1)) x(2:end)], Labels{m}, 'UniformOutput', false);
                % Show legend
                lgd = legend(Labels2{m},"AutoUpdate","off","Location","northeast");
    
                xticklabels(XLab);
                xlabel('Cluster Number')
                ylabel('Number of events')
                xLimits = xlim;
                yLimits = ylim;
                text(0.5*xLimits(2),0.95*yLimits(2),T_Variables{m},"HorizontalAlignment","center","VerticalAlignment","middle","FontWeight","bold")
                text(0.02*xLimits(2),0.95*yLimits(2),[letter,'.'],"HorizontalAlignment","left","VerticalAlignment","middle","FontWeight","bold")
                letter = char(letter+1);
    
                XLimitsBar = xLimits;
    %             ax = gca;
    %             ax.FontSize = 20;
    %             lgd.FontSize = 14;
    
            end
            
            if Ground_Truth_Figs ==1
                AllClasses = unique(MyClasses(:,1));
                AllClasses = string(AllClasses{:,1});
        
                % Get events of each class in each cluster
                for k = 1:numclust
                    ClustIdx = find(clustI==k);
                    MyClass{k} = string(MyClasses{ClustIdx,:});
                    for j = 1:length(AllClasses)
                        NumClass(k,j) = nnz(strcmp(MyClass{k},AllClasses{j}));
                    end
                end
        
                % Create simple figure
                nexttile
                bx = bar(NumClass(sortedIdx,:),'stacked','FaceColor','flat');
                % Change colours to be colorblind friendly
                NewCols = [Colorblind2; 0 0 0];
                for x = 1:length(bx)
                    bx(x).CData = NewCols(x,:);
                end
                xLimits = xlim;
                yLimits = ylim;
                GT_Label = MyClasses.Properties.VariableNames;
                text(0.5*xLimits(2),0.95*yLimits(2),GT_Label,"HorizontalAlignment","center","VerticalAlignment","middle","FontWeight","bold")
                text(0.02*xLimits(2),0.95*yLimits(2),[letter,'.'],"HorizontalAlignment","left","VerticalAlignment","middle","FontWeight","bold")
                ylabel('Number of events')
                xlabel('Cluster number')
                xticklabels(XLab)
                lgd_bx = legend(bx,AllClasses,'Location','northeast','AutoUpdate','Off');
            end
        end
    end

    %% Plot box plots of cluster statistics

        for k = 1:length(numclust)
            % Get cluster indexes for numclust(k)
            clustI = clust(:,k); 
        
            % Initiate figure
            fig = figure('WindowStyle','docked');
            numRows = ceil(NumVars/2);
            TL = tiledlayout(numRows,2,'TileSpacing','compact','Padding','tight');
    
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
    
            ClustAdjusted{k} = clustI;
    
            % Work out data rows of outlier clusters
            IdxOutlier = ClustCounts == 1;
            SingleClusts = find(IdxOutlier);
            SingleClustOrder = sort(SingleClusts,'ascend'); % To make appending below simpler
    
            col = colorTriplets(2:end,:); % To avoid black as the first colour
            if ~isempty(SingleClusts)
                for j = 1:length(SingleClusts)
                    IdxOutlierData(j) = find(clustI == SingleClusts(j));
                
                    % Add black to color triplets for outliers so colours of
                    % violins match dendrogram
                    col = [col(1:SingleClustOrder(j)-1,:); [0,0,0]; col(SingleClustOrder(j):end,:)];
                end
            end
            % Cut colour array to correct size for number of clusters
            col = col(1:numclust(k),:);
    
            ParamNamesUnits = ParamNames;
    
            % Plot box plot for each deformation parameter
            letter = 'a';
            for i = 1:NumVars
                nexttile
                hold on
                
                OrigData = TCropOrig(:,i);
                LogData = TCropLog(:,i);
                Data = TCrop(:,i);
    
                %vp = daviolinplot(Data,'groups',clustI,'violin','half','colors',col,'box',2,'boxcolors','same','outliers',1);
                vp = daboxplot(Data,'groups',clustI,'fill',1,'colors',col,'whiskers',1,'scatter',2,'scattersize',40,'flipcolors',0,'jitter',1,'mean',1,'outliers',1,'outsymbol','k*','boxalpha',0.75);
                box on
                xticklabels(XLab);
                xlabel('Cluster number')
                %set(gca,'XTick',[])
                ylabel('Feature value')
    
                % Delete violins of single-event clusters
                delete(vp.bx(IdxOutlier));
                delete(vp.md(IdxOutlier));
                delete(vp.mn(IdxOutlier));
                delete(vp.ot(IdxOutlier));
                delete(vp.wh(1,IdxOutlier,:));
                delete(vp.sc(IdxOutlier));
                if ~isempty(SingleClusts)
                    scatter(SingleClusts,Data(IdxOutlierData),36,'k','Marker','_','LineWidth',2);
                end
    
                xlim(XLimitsBar);
    
                drawnow
    
                xLimits = xlim;
                yLimits = ylim;
                YPos = yLimits(1) + 0.95*(yLimits(2)-yLimits(1));
                text(0.5*xLimits(2),YPos,ParamNamesUnits{i},"HorizontalAlignment","center","VerticalAlignment","middle","FontWeight","bold")
                text(0.02*xLimits(2),YPos,[letter,'.'],"HorizontalAlignment","left","VerticalAlignment","middle","FontWeight","bold")
                letter = char(letter+1);
    
                % Try to scale non-transformed values
                % Find range of original values

                if AddRight ==1
                    origMin = min(OrigData);
                    origMax = max(OrigData);
                    origTicks = niceTicks(origMin,origMax,5);

                    % Transforming data
                    if Lg ==1 % Optionally take log of absolute values
                        tickPos = abs(origTicks);
                        tickPos = log10(origTicks);
                    elseif Lg==2 % Take log and preserve sign
                        Signs = sign(origTicks);
                        tickPos = abs(origTicks);
                        tickPos = log10(1+origTicks);
                        tickPos = origTicks.*Signs;
                    elseif Lg==3 % Inverse hyperbolic sign
                        tickPos = asinh(origTicks);
                    end

                    % Normalising data
                    if Normalisation ==1
                        % ZScore
                        %tickPos = zscore(tickPos);
                        tickPos = (tickPos-mu)./sigma;
                    elseif Normalisation ==2
                        % Take log
                        tickPos = abs(tickPos);
                        tickPos = log10(tickPos);
                    elseif Normalisation ==3
                        % Scale between -1 and 1
                        MaxVal=max(abs(LogData));
                        tickPos = tickPos./MaxVal;
                    end

                    % Find equivalent feature values for YTicks

                    yyaxis right
                    ylim(yLimits);
                    yticks(tickPos);
                    yticklabels(string(origTicks));
                    ylabel('Original value');
        
                    ax = gca;
                    ax.YAxis(1).Color = [0 0 0];
                    ax.YAxis(2).Color = [0 0 0];
                end
    
    
                % Add legend
                if i==NumVars
                    hold on
                    % Dummy plot
                    for j = 1:(length(ClustCounts)-length(SingleClusts))
                        violgdplot(j) = plot(nan,nan,'Color',colorTriplets(j+1,:),'DisplayName',num2str(j),'LineWidth',3,'LineStyle','-');
                    end
                    violgdplot(end + 1) = scatter(nan,nan,36,'k','Marker','_','DisplayName','Outliers','LineWidth',3);
    
                    lgd = legend(violgdplot,'Location','eastoutside','AutoUpdate','off');
                    title(lgd,'Cluster');
                    box on
                end
            end
        end

        %% Slihouette score plot
        figure('WindowStyle','docked')
        Tl2 = tiledlayout(1,length(numclust));
        for k = 1:length(numclust)
            nexttile
            hold on
            [SVals, SFig] = silhouette(TCrop,clustI,DistMethod)
            %silhouette(TCrop,clust(:,k),DistMethod)
            s2(k) = mean(silhouette(TCrop,clustI,DistMethod));
            s3(k) = median(silhouette(TCrop,clustI,DistMethod));
            view([90 -90])
            set(gca,'YDir','normal')
            box on

            % Get axes
            ax = SFig.Children;
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
            xline(0.5,':')
            SilBar.CData = C;
            SilBar.FaceColor = 'flat';
        end

    %% Do summary figures based on silhouette score
    for k = 1:length(numclust)
        figure('WindowStyle','Docked')
        Tl2 = tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
        for i = 1:numclust
            IdxClust2 = clustI==i;
            for j = 1:NumVars
                nexttile(j)
                hold on
                scat(i) = scatter(SVals(IdxClust2),...
                    TCrop(IdxClust2,j), ...
                    50, ...
                    col(i,:), ...
                    'filled', ...
                    'MarkerFaceAlpha',0.75);
    
                ylabel('Feature Value')
                yLimits = [min(TCrop(:,j)), max(TCrop(:,j))];
                ylim(yLimits);
                
                yyaxis right
                ylim(yLimits);
                yticks(tickPos);
                yticklabels(string(origTicks));
                ylabel('Original Value')
                ax = gca;
                ax.YAxis(1).Color = [0 0 0];
                ax.YAxis(2).Color = [0 0 0];
    
                xlabel('Silhouette score')
                box on
                if i==numclust
                    text(0.5,0.95,ParamNamesUnits{j}, ...
                    'Units','normalized', ...
                    'HorizontalAlignment','center', ...
                    'VerticalAlignment','top', ...
                    'FontWeight','bold');
                end
            end
        end
        % Add legend
        nexttile(1)
        leg_sil = legend(scat, XLab, 'Location', 'northwest');
        title(leg_sil,'Cluster')
        
        %% Save Cluster stats to Excel table
        ClustVals = unique(clustI);
        for i = 1:length(ClustVals)
            for j = 1:width(TCrop)
                IdxThis = clustI==ClustVals(i);
            
                numVals = sum(IdxThis);
                meanVals = mean(TCropOrig(IdxThis,j));
                medVals = median(TCropOrig(IdxThis,j));
                maxVals = max(TCropOrig(IdxThis,j));
                minVals = min(TCropOrig(IdxThis,j));
                PC_25 = prctile(TCropOrig(IdxThis,j),25);
                PC_75 = prctile(TCropOrig(IdxThis,j),75);
                IQR = PC_75-PC_25;
                SValsMed = median(SVals(IdxThis));
                SValsMean = mean(SVals(IdxThis));
                SValsMax = max(SVals(IdxThis));
                SValsMin = min(SVals(IdxThis));
                S_NonZero = sum(SVals(IdxThis)>0);
                S_Zero = sum(SVals(IdxThis)<0);
                
                ClustStatsFull(j,i,:) = [i numVals meanVals medVals maxVals minVals PC_25 PC_75 IQR SValsMed SValsMean SValsMax SValsMin S_NonZero S_Zero];
            end
        end
        
        ClustStatsTab1 = array2table(squeeze(ClustStatsFull(1,:,:)),'VariableNames',...
            {'Cluster no.', 'Num events', 'Mean', 'Median', 'Max', 'Min', '25%', '75%','IQR','Median S', 'Mean S', 'Max S', 'Min S', 'S>0', 'S<0'});
        ClustStatsTab2 = array2table(squeeze(ClustStatsFull(2,:,:)),'VariableNames',...
            {'Cluster no.', 'Num events', 'Mean', 'Median', 'Max', 'Min', '25%', '75%','IQR','Median S', 'Mean S', 'Max S', 'Min S', 'S>0', 'S<0'});
        ClustStatsTab3 = array2table(squeeze(ClustStatsFull(3,:,:)),'VariableNames',...
            {'Cluster no.', 'Num events', 'Mean', 'Median', 'Max', 'Min', '25%', '75%','IQR','Median S', 'Mean S', 'Max S', 'Min S', 'S>0', 'S<0'});
        ClustStatsTab4 = array2table(squeeze(ClustStatsFull(4,:,:)),'VariableNames',...
            {'Cluster no.', 'Num events', 'Mean', 'Median', 'Max', 'Min', '25%', '75%','IQR','Median S', 'Mean S', 'Max S', 'Min S', 'S>0', 'S<0'});
        ClustStatsTab5 = array2table(squeeze(ClustStatsFull(5,:,:)),'VariableNames',...
            {'Cluster no.', 'Num events', 'Mean', 'Median', 'Max', 'Min', '25%', '75%','IQR','Median S', 'Mean S', 'Max S', 'Min S', 'S>0', 'S<0'});
        
        ClusterTabs = {ClustStatsTab1,ClustStatsTab2,ClustStatsTab3,ClustStatsTab4,ClustStatsTab5};
        
        if SaveClustStats ==1
            mkdir('ClusterStats')
            filename = strcat(OutFolder,'/Stats_',RunName,'_',num2str(numclust(k)),'.xlsx');
            filename2 = strcat(OutFolder,'/Clusters_',RunName,'_.xlsx');
            %Write each table to a separate sheet in the Excel sheet
            for i = 1:length(ClusterTabs)
                sheetName = ParamNames{i}; % Change to your desired sheet names
                writetable(ClusterTabs{i}, filename, 'Sheet', sheetName);
            end
            writetable(table(clust(:,k)), filename2, 'Sheet', [num2str(numclust(k)),'_Clusters']);
        end
    end


    %% Helper function
    function ticks = niceTicks(vmin,vmax,nmax)

        span = vmax-vmin;

        step0 = span/(nmax-1);

        mag = 10^floor(log10(step0));

        candidates = [1 2 5 10]*mag;

        [~,idx] = min(abs(candidates-step0));

        step = candidates(idx);

        firstTick = ceil(vmin/step)*step;
        lastTick  = floor(vmax/step)*step;

        ticks = firstTick:step:lastTick;
    end
end