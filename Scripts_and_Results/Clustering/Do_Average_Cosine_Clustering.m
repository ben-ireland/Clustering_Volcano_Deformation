function OutFile = Do_Average_Cosine_Clustering(MyClasses,T,TInfo)
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
    %% Options
    RunName = 'Cosine_Average_Paper';
    
    % Save options
    SaveClusters = 1; % Save clusters to put into coloured dendrogram
    SaveClustStats = 1;
    
    % Figures
    GrpStat = 1; %Run group stats portion of the script?
    DendroLab = 0;  %Display dendrogram labels or not - doesn't work if set to 1
    LongLabels = 1; % Display long labels? (NAME (STARTDATE - ENDDATE))
    AddRight = 1; % Add representative values to Y axis for cluster stats (sometimes doesn't work)
    Comp_MyClasses = 1; % Optionally compare with my classifications for the events

    if Comp_MyClasses ==1
        MyClasses = readtable('/Users/jl20461/Library/CloudStorage/OneDrive-UniversityofBristol/Documents/BristolPhD/Clustering_Tests/June25/Subset_Du_Ra_Ar_As_De_DefCat_All_Erup_InSAR3n.xlsx',Sheet='Classification2_NoDup');
        MyClasses = MyClasses.Classification_Alt2;
    end
    
    numclust = 9;
    DistMethod = 'cosine';
    LinkMethod = 'average';

    %% Input and pre-process deformation catalogue
    T.Properties.VariableNames = {'VolcNum','Name','Duration','Velocity','MaxD','Area','AR','MinDepth','MaxDepth'};
    T{:,end+1} = (T.MinDepth + T.MaxDepth)./2;
    T.Properties.VariableNames{end} = 'MeanDepth';
    T = removevars(T,{'MinDepth','MaxDepth','VolcNum','Name','MaxD'}); % Remove any other variables here
   
    
    %% Creating labels/subsets of string parameters and removing NaNs
    % Remove NaNs
    ParamNames = T.Properties.VariableNames;
    ParamNamesUnits = {'Duration (years)', 'Velocity (cm yr^{-1})','Area (km^2)','Aspect Ratio','Mean Modelled Depth (km)'};
    NumVars = length(ParamNames);
    T = table2array(T);
    NanMask = ~any(isnan(T), 2);
    
    TCrop = T(NanMask,:);
    TInfoCrop = TInfo(NanMask,:);
    
    if LongLabels==1
        labels = append(TInfoCrop.Volcano,' (',TInfoCrop.StartDate,' - ',TInfoCrop.EndDate,')');
    else
        labels = TInfoCrop.Volcano;
    end
    
    % Copy original cropped dataset
    TCropOrig = TCrop;

    % Apply IHS transformation
    TCrop = asinh(TCrop); % inverse hyperbolic sine
    Lg_name = 'log_IHS_';
    
    TCropLog = TCrop;
    TCrop2 = abs(TCropOrig);
    
    % Normalise scores between -1 and 1
    for i = 1:size(TCrop,2)
        MaxVal(i)=max(abs(TCrop(:,i)));
        TCrop(:,i) = TCrop(:,i)./MaxVal(i);
    end

    % Visualise distribution of parameters with histograms
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
    
        if k==1
            title('Duration (years)')
        elseif k==2
            title('Velocity (cm yr^{-1})')
        elseif k==3
            title('Area (km^2)')
        elseif k==4
            title('Aspect ratio')
        elseif k==5
            title('Mean depth (km)')
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

    % Indexes of table columns
    KeptIdx = 1:width(TCrop);
    
    %% Cluster data Hierarchically
    distances=pdist(TCrop,DistMethod);
    links=linkage(distances, LinkMethod);
    
    Folder = pwd;
    if SaveClusters==1
        Z = links;
        X = TCrop;
        mkdir('Paper_Clustering_Results')
        OutFile = [Folder,'/Paper_Clustering_Results/ZX_Data_',RunName,'.mat'];
        save(OutFile,'Z','X');
    end
    
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
        if DendroLab ==1
            [H, ~, P, colors{k},theGroups(:,k)] = dendrogram2(links, 0, 'Labels',labels, 'orientation','left','ColorThreshold',cutoff_for_colours(k));
        elseif DendroLab ==0
            [H, ~, P, colors{k},theGroups(:,k)] = dendrogram2(links, 0, 'Labels',[], 'orientation','left','ColorThreshold',cutoff_for_colours(k));
        end
        hold on

        T = clust(:,k);                 % your cluster() output
        Cmap = colorTriplets(2:end,:);
        NumMultClusts = recolorDendrogramClasses(H, T, labels, Cmap, 1);

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
        ylabel('Volcano')
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

    
    %     savefig(['Dendrogram_',num2str(numclust(k)),Lg_name,'.fig'])
    %              saveas(gcf,['Dendrogram_',num2str(numclust(k)),'_LogNorm.png'],'png')
    %              saveas(gcf,['Dendrogram_',num2str(numclust(k)),'_LogNorm.svg'],'svg')
            
    end
    
   
    
    %% Prep stats for each cluster
    if GrpStat ==1
        
        % Process Info table
        %Group together similar groups for sources
        idxSourceNan = isempty(TInfoCrop.SourceOfDeformation);
        TInfoCrop.SourceOfDeformation(idxSourceNan) = {'Unspecified'};
        SourceStr = {'magma (eruption)','magma (non-eruption)','magma (non-eruption); anthropogenic','magma (non-eruption); faulting','magma (non-eruption); hydrothermal','magma (post-eruption)','faulting; magma (eruption)'};
        SourceStr2 = {'hydrothermal','hydrothermal system'};
        SourceStr3 = {'Unknown','Unspecified',''};
    
        idxSource = ismember(TInfoCrop.SourceOfDeformation,SourceStr);
        idxSource2 = ismember(TInfoCrop.SourceOfDeformation,SourceStr2);
        idxSource3 = ismember(TInfoCrop.SourceOfDeformation,SourceStr3);
        TInfoCrop.SourceOfDeformation(idxSource) = {'magma'};
        TInfoCrop.SourceOfDeformation(idxSource2) = {'hydrothermal'};
        TInfoCrop.SourceOfDeformation(idxSource3) = {'Unspecified'};
    
        % Merge volcano types
        VolcTypeStr = {'Stratovolcano','Stratovolcano(es)'};
        VolcTypeStr2 = {'Caldera','Caldera(s)'};
        VolcTypeStr3 = {'Shield','Shield(s)'};
        idxVT = ismember(TInfoCrop.VolcanoType,VolcTypeStr);
        idxVT2 = ismember(TInfoCrop.VolcanoType,VolcTypeStr2);
        idxVT3 = ismember(TInfoCrop.VolcanoType,VolcTypeStr3);
        TInfoCrop.VolcanoType(idxVT) = {'Stratovolcano(es)'};
        TInfoCrop.VolcanoType(idxVT2) = {'Caldera(s)'};
        TInfoCrop.VolcanoType(idxVT3) = {'Sheild(s)'};
    
        % Merge major rock types
        RockStr = {'No Data (checked)','Not on GVP'};
        MaficStr = {'Basalt / Picro-Basalt','Foidite','Trachybasalt / Tephrite Basanite'};
        IntStr = {'Andesite / Basaltic Andesite','Dacite'};
        FelsicStr = {'Phonolite','Rhyolite','Trachyte / Trachydacite'};
    
        idxRock = ismember(TInfoCrop.Majorrocktype,RockStr);
        idxMaf = ismember(TInfoCrop.Majorrocktype,MaficStr);
        idxInt = ismember(TInfoCrop.Majorrocktype,IntStr);
        idxFel = ismember(TInfoCrop.Majorrocktype,FelsicStr);
        TInfoCrop.Majorrocktype(idxRock) = {'Unspecified'};
        TInfoCrop.Majorrocktype(idxMaf) = {'Mafic'};
        TInfoCrop.Majorrocktype(idxInt) = {'Intermediate'};
        TInfoCrop.Majorrocktype(idxFel) = {'Felsic'};
    
        % Simplify tectonic setting
        RiftStr = {'Rift zone'};
        IntrStr = {'Intraplate'};
        SubdStr = {'Subduction'};
    
        idxRift = contains(TInfoCrop.Tectonicsetting,RiftStr);
        idxIntr = contains(TInfoCrop.Tectonicsetting,IntrStr);
        idxSubd = contains(TInfoCrop.Tectonicsetting,SubdStr);
    
        TInfoCrop.Tectonicsetting(idxRift) = {'Rift zone'};
        TInfoCrop.Tectonicsetting(idxIntr) = {'Intraplate'};
        TInfoCrop.Tectonicsetting(idxSubd) = {'Subduction zone'};
    
        % Linked eruption
        TInfoCrop.Eruption_ = string(TInfoCrop.Eruption_);
        TInfoCrop.Eruption_(TInfoCrop.Eruption_ == "1") = "No linked eruption";
        TInfoCrop.Eruption_(TInfoCrop.Eruption_ == "2") = "Linked eruption";
    
        % Rename variables and remove ones we're not using
        T = TInfoCrop;
        T = removevars(T,{'VolcanoNumber_GVP_','Volcano','Location','TypeOfRate_mean_Max_Estimated_','Direction_original_'});
        T_Variables = {'Linked Eruption','Volcano Type','Tectonic Setting','Major Rock Type','Deformation direction','Deformation Source','Start Date','End date'};
        T.Properties.VariableNames = T_Variables;
    
        nCols = width(T)-2;     
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
    
        varNames = T_Variables(1:end-2);
        allCats = unique(allCats); 
        numCats = numel(allCats);
    
        C = zeros(numCats, nCols);   % rows = categories, columns = table columns
        for col = 1:nCols
            colData = T.(varNames{col});
            for k = 1:numCats
                C(k, col) = sum(strcmp(colData, allCats{k}));
            end
        end
    
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
        ylim([0 height(TInfoCrop)]);
        yLims = ylim;
        xticks(1:nCols);
        xticklabels(varNames);
    
        yyaxis right
        ylim(yLims)
        TickVals = [0 height(TInfoCrop)./4 height(TInfoCrop)./2 3*(height(TInfoCrop)./4) height(TInfoCrop)];
        yticks(TickVals)
        yticklabels({'0','25','50','75','100'})
        ylabel('Proportion of signals (%)')
        ax = gca;
        ax.YAxis(1).Color = [0 0 0];
        ax.YAxis(2).Color = [0 0 0];
        box on
    
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
    
    
        % Volcanolgical characteristics comparison
        for k = 1:length(numclust)
            clustI = clust(:,k); 
    
            fig = figure('WindowStyle','docked');
            TL = tiledlayout(4,2,'TileSpacing','compact','Padding','tight');
            if length(numclust)>1
                txt2 = title(TL,"Clusters: " + numclust(k));
                txt2.FontSize = 16;
            end
            letter = 'f'; % For subplot numbers
            nexttile % Empty plot
            set(gca,'XTick',[])
            set(gca,'YTick',[])
            ax = gca;
            ax.XAxis.Color = [1 1 1];
            ax.YAxis.Color = [1 1 1];
            box off
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
                %xlabel('Cluster Number')
                if m==2 || m==4 || m==6
                    ylabel('Number of events')
                end
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
            
            AllClasses = unique(MyClasses);
            % Get events of each class in each cluster
            for k = 1:numclust
                ClustIdx = find(clustI==k);
                MyClass{k} = MyClasses(ClustIdx);
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
            text(0.5*xLimits(2),0.95*yLimits(2),'Inferred process',"HorizontalAlignment","center","VerticalAlignment","middle","FontWeight","bold")
            text(0.02*xLimits(2),0.95*yLimits(2),[letter,'.'],"HorizontalAlignment","left","VerticalAlignment","middle","FontWeight","bold")
            %ylabel('Number of events')
            %xlabel('Cluster number')
            xticklabels(XLab)
            lgd_bx = legend(bx,AllClasses,'Location','northeast','AutoUpdate','Off');
        end
    
            
    %        savefig(['Clust_Stats_',num2str(numclust(k)),Lg_name,Norm_name,'.fig'])
    %               saveas(gcf,['Clust_Stats_',num2str(numclust(k)),'_LogNorm.png'],'png')
    %               saveas(gcf,['Clust_Stats_',num2str(numclust(k)),'_LogNorm.svg'],'svg')
        end
    
        % Violin plots of cluster stats
        for k = 1:length(numclust)
            % Get cluster indexes for numclust(k)
            clustI = clust(:,k); 
        
            % Initiate figure
            fig = figure('WindowStyle','docked');
            TL = tiledlayout(4,2,'TileSpacing','compact','Padding','tight');
    %         if length(numclust)>1
    %             txt2 = title(TL,"Clusters: " + numclust(k));
    %             txt2.FontSize = 16;
    %         end
            % Empty plots to make same size as volcano characteristics plots
            nexttile % Empty plot
            set(gca,'XTick',[])
            set(gca,'YTick',[])
            ax = gca;
            ax.XAxis.Color = [1 1 1];
            ax.YAxis.Color = [1 1 1];
            box off
    
            nexttile % Empty plot
            set(gca,'XTick',[])
            set(gca,'YTick',[])
            ax = gca;
            ax.XAxis.Color = [1 1 1];
            ax.YAxis.Color = [1 1 1];
            box off
    
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
    
            ParamNamesUnits = {'Duration (years)', 'Velocity (cm yr^{-1})','Area (km^2)','Aspect Ratio','Mean Modelled Depth (km)'};
    
            % Plot violin for each deformation parameter
            letter = 'a';
            for i = 1:NumVars
                nexttile
                hold on
                
                OrigData = TCropOrig(:,i);
                Data = TCrop(:,i);
    
                %vp = daviolinplot(Data,'groups',clustI,'violin','half','colors',col,'box',2,'boxcolors','same','outliers',1);
                vp = daboxplot(Data,'groups',clustI,'fill',1,'colors',col,'whiskers',1,'scatter',2,'scattersize',40,'flipcolors',0,'jitter',1,'mean',1,'outliers',1,'outsymbol','k*','boxalpha',0.75);
                box on
                xticklabels(XLab);
                if i>3
                    xlabel('Cluster number')
                end
                %set(gca,'XTick',[])
                if i==1 || i==3 || i==5
                    ylabel('Feature value')
                end
    
                % Delete violins of single-event clusters
    %             delete(vp.ds(IdxOutlier));
    %             delete(vp.bx(IdxOutlier));
    %             delete(vp.md(IdxOutlier));
    %             delete(vp.wh(IdxOutlier));
    %             delete(vp.ot(IdxOutlier));
    %             delete(vp.ds(IdxOutlier));
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
                RngData = max(OrigData) - min(OrigData);
                numVals = length(KeptIdx);
    %             DataForVals = [prctile(OrigData,2.5), prctile(OrigData,25), prctile(OrigData,50), prctile(OrigData,75),prctile(OrigData,97.5)];
    %             [GC, GR] = groupcounts(round(OrigData./10));   
    %             DataForVals = round(DataForVals,0); % Values for YTickLabels
    %             DataForVals = unique(DataForVals);
                
                if AddRight ==1
                    if KeptIdx(i)==1
                        DataForVals = [0.5 1 2 5 10 20];
                    elseif KeptIdx(i)==2
                        DataForVals = [-1000 -100 -10 -1 1 10 100 1000];
                    elseif KeptIdx(i)==3
                        DataForVals = [1 5 10 100 1000];
                    elseif KeptIdx(i)==4
                        DataForVals = [1 1.5 2 3 4];
                    elseif KeptIdx(i)==5
                        DataForVals = [0.5 1 2 5 10 20];
                    end
        
                    % Find equivalent feature values for YTicks
                    OrigSort = sort(OrigData);
                    DataSort = sort(Data);
                    Idxs = knnsearch(OrigSort,DataForVals');
                    TickLocations = DataSort(Idxs);
                    TickLocations = unique(TickLocations,'stable');
                    allTickLocs{i} = TickLocations;
                    yyaxis right
                    ylim(yLimits);
                    yticks(TickLocations);
                    yticklabels(string(DataForVals));
                    if i==2 || i==4
                        ylabel('Original value');
                    end
        
                    ax = gca;
                    ax.YAxis(1).Color = [0 0 0];
                    ax.YAxis(2).Color = [0 0 0];
                end
    
                % Get stats of clusters
    
    
                % Add legend
                if i==5
                    nexttile
                    hold on
                    % Dummy plot
                    for j = 1:(length(ClustCounts)-length(SingleClusts))
                        violgdplot(j) = plot(nan,nan,'Color',colorTriplets(j+1,:),'DisplayName',num2str(j),'LineWidth',3);
                    end
                    violgdplot(end + 1) = scatter(nan,nan,36,'k','Marker','_','DisplayName','Outliers','LineWidth',3);
    
                    lgd = legend(violgdplot,'Location','west');
                    title(lgd,'Cluster');
                    set(gca,'XTick',[])
                    set(gca,'YTick',[])
                    ax = gca;
                    ax.XAxis.Color = [1 1 1];
                    ax.YAxis.Color = [1 1 1];
                    box off
                end
            end
        end
  
    
    %% Validate clusters
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
            if j==1
                DataForVals = [0.5 1 2 5 10 20];
            elseif j==2
                DataForVals = [-1000 -100 -10 -1 1 10 100 1000];
            elseif j==3
                DataForVals = [1 5 10 100 1000];
            elseif j==4
                DataForVals = [1 1.5 2 3 4];
            elseif j==5
                DataForVals = [0.5 1 2 5 10 20];
            end
            ylim(yLimits);
            yticks(allTickLocs{j});
            yticklabels(string(DataForVals));
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
    
    % Cluster stats
    ClustVals = unique(clustI);
    for k = 1:length(ClustVals)
        for j = 1:width(TCrop)
            IdxThis = clustI==ClustVals(k);
        
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
            
            ClustStatsFull(j,k,:) = [k numVals meanVals medVals maxVals minVals PC_25 PC_75 IQR SValsMed SValsMean SValsMax SValsMin S_NonZero S_Zero];
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
        filename = strcat(Folder,'/ClusterStats/ClusterStats_',RunName,'.xlsx');
        %Write each table to a separate sheet in the Excel file
        for i = 1:length(ClusterTabs)
            sheetName = ParamNames{i}; % Change to your desired sheet names
            writetable(ClusterTabs{i}, filename, 'Sheet', sheetName);
        end
    end
end

