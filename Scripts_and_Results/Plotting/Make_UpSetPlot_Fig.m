function Make_UpSetPlot_Fig(T,T2,T3,TInfo,GVP,GVP2,compareInterp)
%% Ben Ireland, June 2026, Cluster multi-variate datasets
% School of Earth Sciences, University of Bristol
%
% Script to reproduce Upset plot figure
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
        126 041 084;...
        0 0 0;...
        255 255 255;
        150 150 150;
        50 50 50];
    Colorblind2 = Colorblind2./255;
    
    % Process Info table
    T3.Properties.VariableNames = {'VolcNum','Name','Duration','Velocity','MaxD','Area','AR','MinDepth','MaxDepth'};
    T3{:,end+1} = (T3.MinDepth + T3.MaxDepth)./2;
    T3.Properties.VariableNames{end} = 'MeanDepth';
    T3 = removevars(T3,{'MinDepth','MaxDepth','VolcNum','Name','MaxD'}); % Remove any other variables here
    
    % Remove NaNs
    ParamNames = T3.Properties.VariableNames;
    NumVars = length(ParamNames);
    T3 = table2array(T3);
    NanMask = ~any(isnan(T3), 2);
    TInfoCrop = TInfo(NanMask,:);
    
    %% Process info table
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
    T4 = TInfoCrop;
    T4 = removevars(T4,{'VolcanoNumber_GVP_','Volcano','TypeOfRate_mean_Max_Estimated_','Direction_original_'});
    T_Variables = {'Linked Eruption','Location','Volcano Type','Tectonic Setting','Major Rock Type','Deformation direction','Deformation Source','Start Date','End date'};
    T_Variables = {T_Variables{2}, T_Variables{1}, T_Variables{3:end}};
    T4 = T4(:, [2 1 3:end]);
    T4.Properties.VariableNames = T_Variables;
    nCols = width(T4)-2;
    
    %% Define column sets
    columnSets = {6,11,12,14,27,28,[6 27],[6 11 12],[27 6 11 12],[28 6 11 12],[28 6 11 12 14],[27 6 11 12 14]};
    % 3/27 = rate (absolute/signed); 6 = duration; 11 = area; 12 = aspect ratio; 13:14 = source
    % depth estimate; 15/28 = max. displacement (absolute/signed)
    
    % Process numeric data only
    raw = table2cell(T);
    [numRows, numCols] = size(raw);
    num = nan(numRows, numCols);  % Pre-fill with NaNs
    
    if compareInterp ==1
        raw2 = table2cell(T2); % Convert the table to a cell array
        num2 = num;
    end
    
    for i = 1:numRows
        for j = 1:numCols
            if isnumeric(raw{i,j})
                num(i,j) = raw{i,j};     % Copy numeric values
            elseif islogical(raw{i,j})
                num(i,j) = double(raw{i,j});  % Optionally convert logicals
            end
        end
    end
    
    if compareInterp==1
        for i = 1:numRows
            for j = 1:numCols
                if isnumeric(raw2{i,j})
                    num2(i,j) = raw2{i,j};     % Copy numeric values
                elseif islogical(raw2{i,j})
                    num2(i,j) = double(raw2{i,j});  % Optionally convert logicals
                end
            end
        end
    end
    
    %% Counts for bar chart
    numSets = length(columnSets);
    counts = zeros(numSets, 1);
    
    for i = 1:numSets
        cols = columnSets{i};
        dataSubset = num(:, cols);
        completeRows = all(~isnan(dataSubset), 2); % complete rows = no missing
        counts(i) = sum(completeRows);
    
        % unqiue volcanoes
        NumVolcs(i) = length(unique(num(completeRows,1)));
    end
    
    if compareInterp ==1
        for i = 1:numSets
            cols = columnSets{i};
            dataSubset = num2(:, cols);
            completeRowsInt = all(~isnan(dataSubset), 2); % complete rows = no missing
            counts2(i) = sum(completeRowsInt);
    
            % unqiue volcanoes
            NumVolcsInt(i) = length(unique(num(completeRowsInt,1)));
        end
    end
    
    %% Make map
    % Find subsets (with and without interpolation) with no NaN rows
    dataSubsetInt = num2(:,columnSets{end});
    dataSubset = num(:,columnSets{end});
    completeRows = all(~isnan(dataSubset), 2);
    completeRowsInt = all(~isnan(dataSubsetInt), 2);
    
    % Load GVP holocene and pleistocene databases
    GVP_All = [GVP; GVP2];
    GVP_Nums = GVP_All.VolcanoNumber;
    
    % Collect all GVP nums from table
    AllNums = T.VolcanoNumber_GVP_;
    Complete_Nums = T.VolcanoNumber_GVP_(completeRows);
    Complete_NumsInt = T.VolcanoNumber_GVP_(completeRowsInt);
    
    % Find idx of complete rows within the table
    idxComplete = find(completeRows);
    idxCompleteInt = find(completeRowsInt);
    
    % Check if any volcano nums are not found in the GVP database
    IdxGVP = ismember(Complete_Nums,GVP_Nums);
    IdxNo_complete = find(IdxGVP==0);
    IdxNo_global = idxComplete(IdxNo_complete);
    
    IdxGVPInt = ismember(Complete_NumsInt,GVP_Nums);
    IdxNo_completeInt = find(IdxGVPInt==0);
    IdxNo_globalInt = idxCompleteInt(IdxNo_completeInt);
    
    % Find which Idx of GVP database correspond to volcano nums from tables
    IdxGVP2 = ismember(GVP_Nums,Complete_Nums);
    IdxGVP2Int = ismember(GVP_Nums,Complete_NumsInt);
    IdxGVPAll = ismember(GVP_Nums,AllNums);
    
    % Find corresponding rows
    GVPIdx = find(IdxGVP2);
    GVPIdxInt = find(IdxGVP2Int);
    GVPIdxAll = find(IdxGVPAll);
    
    % Find only rows that were completed with interpolation
    DiffIdx = find(~ismember(GVPIdx,GVPIdxInt));
    
    % Extract the lat-lon values for relevant rows from GVP database
    LatLon = [GVP_All.Latitude(GVPIdx(DiffIdx)), GVP_All.Longitude(GVPIdx(DiffIdx))];
    LatLonInt = [GVP_All.Latitude(GVPIdxInt), GVP_All.Longitude(GVPIdxInt)];
    LatLonInt(end+1,:) = [19.742 68.691]; % Add for 1 volcano missing from GVP
    LatLonAll = [GVP_All.Latitude, GVP_All.Longitude];
    LatLonAllDef = [GVP_All.Latitude(GVPIdxAll), GVP_All.Longitude(GVPIdxAll)];
    
    %% start figure
    figure('WindowStyle','docked');
    ExtraRows = 3;
    t = tiledlayout(3,4,"TileSpacing","tight","Padding","tight","TileIndexing","columnmajor");
    count = 'a';
    
    %% Dot matrix for column sets
    % Matrix of presence per column in each set (for dot display)
    allCols = unique([columnSets{:}]);
    ParamLabels = T.Properties.VariableNames(allCols);
    ParamLabels = {'Duration','Signal area','Aspect ratio','Souce depth','Deformation rate','Max. Displacement'};
    numAllCols = length(allCols);
    presenceMatrix = zeros(numSets, numAllCols);
    
    for i = 1:numSets
        presenceMatrix(i, ismember(allCols, columnSets{i})) = 1;
    end
    
    %% Make figure
    
    % Plot bar chart of counts
    nexttile([1 4])
    yyaxis left
    b = bar(counts, 'FaceColor', [240 228 66]./255);
    if compareInterp ==1
        hold on;
        b2 = bar(counts2, 'FaceColor', [0 158 115]./255);
        legend([b b2],{'With interpolation','No interpolation'},'AutoUpdate','off');
        hold off;
    end
    xlim([0 length(columnSets)+1])
    %xticklabels(compose('Set %d', 1:numSets));
    xticklabels('');
    ylabel('Number of Signals (bars)');
    
    yyaxis right
    hold on
    for i = 1:numSets
        plot(i,NumVolcs(i),'-^','MarkerFaceColor',[240 228 66]./255,'MarkerEdgeColor','k','MarkerSize',20)
        plot(i,NumVolcsInt(i),'-^','MarkerFaceColor',[0 158 115]./255,'MarkerEdgeColor','k','MarkerSize',20)
    end
    ylabel('Number of volcanoes (points)')
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';
    
    % Plot subplot letter
    drawnow
    xLimits = xlim;
    yLimits = ylim;
    text(0.02*xLimits(2), ...
     0.98*yLimits(2), ...
     strcat(count,'.'), 'FontWeight', 'bold', 'FontSize', 12, ...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    
    % Plot dot matrix
    nexttile([1 4])
    hold on
    for i = 1:length(columnSets) % Num. unique parameters
        for j = 1:length(ParamLabels) % Num. column sets to compare
            if presenceMatrix(i,j)==1
                if i ~= length(columnSets)
                    % If included in this set, plot it filled
                    plot(i,j,LineStyle="none",Marker="o",MarkerFaceColor='k',MarkerEdgeColor='k',MarkerSize=15);
                else
                    plot(i,j,LineStyle="none",Marker="o",MarkerFaceColor='r',MarkerEdgeColor='r',MarkerSize=15);
                end
            else
                % If not, plot it unfilled
                plot(i,j,LineStyle="none",Marker="o",MarkerFaceColor='none',MarkerEdgeColor='k',MarkerSize=15);
            end
        end
        % Plot lines between params in each column set
        if length(columnSets{i})>1
            if i ~= length(columnSets)
                plot([i,i],[find((presenceMatrix(i,:)==1),1,'first'),find((presenceMatrix(i,:)==1),1,'last')],LineStyle="-",Color='k',Marker='none');
            else
                plot([i,i],[find((presenceMatrix(i,:)==1),1,'first'),find((presenceMatrix(i,:)==1),1,'last')],LineStyle="-",Color='r',Marker='none');
            end
        end
    end
    % Plot lines between 
    % Format axes
    box on
    xlim([0 length(columnSets)+1]);
    ylim([0 numAllCols+1]);
    yticks(1:length(ParamLabels));
    yticklabels(ParamLabels)
    xticks(1:length(columnSets))
    %xticklabels(compose('Set %d', 1:numSets));
    xticklabels('');
    %xlabel('Parameter sets')
    
    % % Plot subplot letter
    % drawnow
    % xLimits = xlim;
    % yLimits = ylim;
    % text(0.02*xLimits(2), ...
    %  0.98*yLimits(2), ...
    %  strcat(count,'.'), 'FontWeight', 'bold', 'FontSize', 12, ...
    %  'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    
    hold off
    
    %% Final plot   
    allCats = {};
    for c = 1:nCols
        % Extract the column as an array
        colData = T4{:, c};
    
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
        colData = T4.(varNames{col});
        for k = 1:numCats
            C(k, col) = sum(strcmp(colData, allCats{k}));
        end
    end
    
    %figure('WindowStyle','docked');
    %tlStats = tiledlayout(nCols+3,nCols,"Padding","tight","TileSpacing","tight",'TileIndexing','rowmajor');
    
    Extra2 = 3;
    nexttile([1 4])
    
    count = char(count+1);
    
    % Main axes
    ax1 = gca;
    hold(ax1,'on')
    
    for k = 1:nCols
        
        SB_Data{k} = C(C(:,k)>0,k);
    
        % Horizontal stacked bars
        h{k} = barh(k, SB_Data{k}, 'stacked');
    
        for a = 1:length(SB_Data{k})
            h{k}(a).FaceColor = Colorblind2(a,:);
        end
    end
    
    % ---- Primary axis (counts) ----
    xlabel(ax1,'Number of events')
    ylim(ax1,[0.5 nCols+0.5])
    xlim(ax1,[0 height(TInfoCrop)])
    
    yticks(ax1,1:nCols)
    yticklabels(ax1,varNames)
    
    xLims = xlim(ax1);
    
    drawnow
    
    % Panel label
    text( ...
        0.02*xLims(2), ...
        0.98*nCols, ...
        strcat(count,'.'), ...
        'FontWeight','bold', ...
        'FontSize',12, ...
        'VerticalAlignment','top', ...
        'HorizontalAlignment','left', ...
        'BackgroundColor','w' );
    
    box(ax1,'on')
    
    % % =========================================================
    % % Secondary x-axis (percentages)
    % % =========================================================
    % 
    % ax2 = axes( ...
    %     'Position', ax1.Position, ...
    %     'Color', 'none', ...
    %     'XAxisLocation', 'top', ...
    %     'YAxisLocation', 'right', ...
    %     'YTick', [], ...
    %     'Box', 'off');
    % 
    % % Match limits
    % ax2.XLim = ax1.XLim;
    % 
    % % Percentage ticks
    % TickVals = [0 ...
    %             height(TInfoCrop)/4 ...
    %             height(TInfoCrop)/2 ...
    %             3*height(TInfoCrop)/4 ...
    %             height(TInfoCrop)];
    % 
    % ax2.XTick = TickVals;
    % ax2.XTickLabel = {'0','25','50','75','100'};
    % 
    % xlabel(ax2,'Proportion of signals (%)')
    % 
    % % Make axes black
    % ax1.XColor = [0 0 0];
    % ax1.YColor = [0 0 0];
    % ax2.XColor = [0 0 0];
    
    figure('WindowStyle','docked')
    tl2 = tiledlayout(Extra2,nCols,'TileSpacing','none','Padding','tight');
    for colnum = 1:nCols
        nexttile([Extra2 1])
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
end
