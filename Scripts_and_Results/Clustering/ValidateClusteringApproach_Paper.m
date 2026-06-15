function OutFolder = ValidateClusteringApproach_Paper(T,TInfo)
%% Ben Ireland, June 2026, Cluster multi-variate datasets
% School of Earth Sciences, University of Bristol
%
% Script to validate different clustering approaches.
%
% Methods in these scripts are described in:
% "Clustering global volcano deformation
% datasets: insights and limitations for analogue signals"
% Submitted to JVGR
%
% Copyright: Ben Ireland, 2026
% V1.0    

%% Colourmaps
    load colorblind_colormap.mat
    load UniqueColors.mat
    IBM_ColorBlind = ...
        [255 176 0;...
        254 97 0;...
        220 38 127;...
        120 94 240;...
        100 143 255];
    IBM_ColorBlind = IBM_ColorBlind./255;
    
    %% Options
    % 1= log values before clustering, 0= don't log values before clustering
    Lg = 3;
    
    % 1 = normalise variables, 2 = normalise by taking the log of the variables, 0 = don't normalise variables
    % 3 = scaled/absolute value
    % value
    Normalisation = 3;
    % 1 = Cluster with PCA; 0 = Cluster in non-PCA space
    PCA = 0;
    NoOutliers=0;
    MinClust = 2;
    MaxClust = 25;
    ClustGap = 1;
    ClustNums = MinClust:ClustGap:MaxClust;
    clustNames = string(MinClust:ClustGap:MaxClust);
    %% Input and data pre-processing for Clustering
    %Ben Ireland, April 2023, cluster deformation data using hierarchicial
    %clustering
    
    T.Properties.VariableNames = {'VolcNum','Name','Duration','Velocity','MaxD','Area','AR','MinDepth','MaxDepth'};
    T{:,end+1} = (T.MinDepth + T.MaxDepth)./2;
    T.Properties.VariableNames{end} = 'MeanDepth';
    T = removevars(T,{'MinDepth','MaxDepth','VolcNum','Name','MaxD'}); % Remove any other variables here
    
    %% Creating labels/subsets of string parameters and removing NaNs
    % Remove NaNs
    ParamNames = T.Properties.VariableNames;
    T = table2array(T);
    NanMask = ~any(isnan(T), 2);
    
    if NoOutliers==1
        OutMask = isoutlier(T,"percentiles",[2.5 97.5],1);
        OutMaskRow = sum(OutMask,2)>0;
    
        TCrop = T(NanMask & ~OutMaskRow,:);
        TInfoCrop = TInfo(NanMask & ~OutMaskRow,:);
    else
        TCrop = T(NanMask,:);
        TInfoCrop = TInfo(NanMask,:);
    end
    
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
    
    if Normalisation ==1
        % ZScore
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
    
    % % Do PCA on the normalised values
    % [coeff,score,latent,tsquared,explained,mu] = pca(TCrop);
    % 
    % % Visualise PCA results
    % figure('WindowStyle','docked')
    % t = tiledlayout(7,6,"TileSpacing","none","Padding","compact");
    % title(t,'PCA coefficients and explained variance')
    % nexttile([6 6])
    % h = heatmap(cellstr(string(1:size(TCrop,2))),ParamNames,abs(coeff),"ColorScaling","scaledrows");
    % ylabel('Deformation parameter')
    % h.XDisplayLabels = repmat({''}, size(h.XDisplayData));
    % hs = struct(h);
    % ylabel(hs.Colorbar, 'PC coefficient (absolute)');
    % colormap(h,"pink")
    % 
    % nexttile([1 6])
    % h2 = heatmap(cellstr(string(1:size(TCrop,2))),'% Variance explained',explained');
    % colorbar off
    % xlabel('Principal component')
    % colormap(h2,"summer")
    % 
    % 
    % figure('WindowStyle','docked')
    % Group = TInfoCrop.Direction_modified_(1:end);
    % % Create full gplotmatrix (all subplots)
    % [h,ax,bigax] = gplotmatrix(TCrop,[],Group,[],[],10,[],'grpbars',ParamNames);
    
    % % Number of variables
    % nVar = size(TCrop,2);
    % 
    % % Hide the upper triangle
    % for i = 1:nVar
    %     for j = 1:nVar
    %         if j > i  % upper triangle
    %             set(ax(i,j),'Visible','off');
    %             delete(get(ax(i,j),'Children')); % remove plot contents
    %         end
    %     end
    % end
    % title(bigax,'Trade-off plots')
    % 
    % % Optional: clean up spacing
    % set(bigax,'Visible','off');
    % set(ax(:),'Box','on');
    
    %% Iteratively testing different distance and linkage methods
    
    distance_method = {'cosine','euclidean','cityblock'};
    metric_labels = {'Cosine', 'Euclidean', 'Cityblock'};
    linkage_method = {'single','average','complete','ward'};
    linkageNames = {'Single', 'Average', 'Complete', 'Ward'};
    
    nRun = 0;
    figure('WindowStyle','docked')
    t = tiledlayout(length(linkage_method),length(distance_method),'TileIndexing','columnmajor','TileSpacing','tight','Padding','tight');
    %title(t,'Dendrograms')
    
    count = 'a';
    for a=1:length(distance_method)
        Distances_Methods(a,:) = pdist(TCrop,distance_method{a});
        
        for n=1:length(linkage_method)
            Links_Methods(:,:,n) = linkage(Distances_Methods(a,:),linkage_method{n});
            [c_M(a,n), cDist(a,n,:)] = cophenet(Links_Methods(:,:,n),Distances_Methods(a,:));
            
            for d=1:1:10
                b = d;
                e(b) = d;
                I(:,:,a,n,b) = inconsistent(Links_Methods(:,:,n),d);
                
                % a = distance method, n = linkage, b = inconsistency depth
                maxi_M(a,n,b)=max(I(:,4,a,n,b));
                meani_M(a,n,b)=mean(I(:,4,a,n,b));
                vari_M(a,n,b)=var(I(:,4,a,n,b));
                stdi_M(a,n,b)=std(I(:,4,a,n,b));
            end
            msg = lastwarn
            if ~contains(msg,'non-Euclidean')
                nRun = nRun + 1;
                
                meanInconsist(nRun,:) = meani_M(a,n,:);
                CopheneticCoeff(nRun) = c_M(a,n);
                CophenetDists(nRun,:) = cDist(a,n,:);
                PairWiseDists(nRun,:) = Distances_Methods(a,:);
    
                if nRun>=4
                    nexttile(nRun+1)
                else
                    nexttile(nRun)
                end
    
                Z = Links_Methods(:,:,n);
                X = TCrop;
                mkdir('Cluster_Validation_Results')
                OutFolder = [pwd,'/Cluster_Validation_Results'];
                save([pwd,'/Cluster_Validation_Results/ZX_Data_',metric_labels{a},'_',linkageNames{n},'.mat'],'Z','X');
    
                % Plot dendorgram for each linkage/distance parameter
                [H, T, P] = dendrogram2(Links_Methods(:,:,n), 0, 'orientation','left');
                hold on
                box on; set(gca,'layer','top');
                RunName{nRun} = strcat(distance_method{a},',',linkage_method{n});
    
                set(H, 'Color', IBM_ColorBlind(n,:)); % Colour dendrogram according to linkage method
                set(H, 'LineWidth', 0.75);
    
                set(gca,'ytick',[]) % Remove data labels
                if nRun==3 || nRun == 7 || nRun == 10 % Add xlabel for bottom row
                    xlabel('Distance')
                end
               
                % Add subplot letter
                xLimits = xlim;
                xlim([0 xLimits(2)*1.2]);
                xLimits = [0 xLimits(2)*1.2]; % Adjust xlims so letter doesn't appear on top
                yLimits = ylim;
                text(0.98*xLimits(2), ...
                 0.98*yLimits(2), ...
                 strcat(count,'.'), 'FontWeight', 'bold', 'FontSize', 12, ...
                 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    
                % Add distance parameter letter
    %             text(0.98*xLimits(2), ...
    %              0.02*yLimits(2), ...
    %              metric_labels{a}, 'FontWeight', 'normal', 'FontSize', 12, ...
    %              'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    
                % Add legend
                if nRun==1
                    hold on
                    for i = 1:length(linkageNames)
                        p2(i) = plot(nan,nan,'LineStyle','-','LineWidth',3,'Marker','none','Color',IBM_ColorBlind(i,:),'DisplayName',linkageNames{i});
                    end
                    lgd3 = legend(p2,'Location','west','AutoUpdate','off');
                    title(lgd3,"Linkage criteria")
                end
    
                % Plot cut-off for varying numbers of clusters as the
                % background of the plot
                Z = Links_Methods(:,:,n);
                hMin = 1e-9;
                hMax = max(Z(:,3));
                nPoints = 1000;  % resolution of the background
                cutHeights = linspace(hMin, hMax, nPoints);
                
                % For each cut height, compute number of clusters
                numClusters = arrayfun(@(h) max(cluster(Z, 'Cutoff', h, 'Criterion', 'distance')), cutHeights);
    
                % Limit max number of clusters for improved visualisation
                numClusters(numClusters>25)=25;
                
                % Create grayscale colormap (avoid pure white)
                % darker for small k, lighter for large k
                grayVals = rescale(numClusters, 0.1, 0.8);
                grayImg = repmat(grayVals, [100 1]);  
                
                % Plot background image
                Img = imagesc([hMin hMax], [0 size(TCrop,1)], grayImg);  % span the same x-range as dendrogram heights
                %Img.AlphaData = 0.5;
                colormap(flipud(gray(length(unique(numClusters)))));
                uistack(H, 'top');
    
                if nRun==10 % Add colourbar to final plot
                    cb = colorbar;
                    cb.Label.String = 'Number of clusters';
                    TickVals = [1,5:5:25];
                    cb.Ticks = rescale(TickVals, 0.1, 0.8);
                    Lab = string(TickVals);
                    Lab(end) = ">25";
                    cb.TickLabels = Lab;
                    cb.Location = "eastoutside";
                end
    
                count = char(count+1);
    
                for k = 1:length(ClustNums)
                    clust = cluster(Links_Methods(:,:,n),'maxclust',ClustNums(k));
    
                    siIndex(k,nRun) = mean(silhouette(TCrop,clust,distance_method{a}));
                    dbIndex(k,nRun) = evalclusters(TCrop, clust, 'DaviesBouldin').CriterionValues;
                    chIndex(k,nRun) = evalclusters(TCrop, clust, 'CalinskiHarabasz').CriterionValues;
                    eva = evalclusters(TCrop, 'linkage', 'gap', 'KList', ClustNums(k));
                    bestK = eva.OptimalK;
                    gapValue(k,nRun) = eva.CriterionValues;
                end
            end
            warning('No warning');
        end
    end
    
    %% Visualise results and validation stats
    
    % if length(ClustNums)<20 % For visualtion clarity and speed
    %     % Plot cluster validation stats
    %     figure('WindowStyle','docked')
    %     t = tiledlayout(1,4,"TileSpacing","loose","Padding","compact");
    %     title(t,'Clustering metrics');
    % 
    %     nexttile
    %     h1 = heatmap(RunName,clustNames,siIndex,"ColorScaling","scaledcolumns");
    %     xlabel('Clustering setup')
    %     ylabel('No. of Clusters')
    %     h1.title('Silhouette coefficient')
    % 
    %     nexttile
    %     h2 = heatmap(RunName,clustNames,dbIndex,"ColorScaling","scaledcolumns");
    %     xlabel('Clustering setup')
    %     h2.title('Davies Bouldin index')
    % 
    %     nexttile
    %     h3 = heatmap(RunName,clustNames,chIndex,"ColorScaling","scaledcolumns");
    %     xlabel('Clustering setup')
    %     h3.title('Calinski Harabasz index')
    % 
    %     nexttile
    %     gapDiff = diff(gapValue);
    %     gapDiff = [zeros(1,length(RunName));gapDiff];
    %     h4 = heatmap(RunName,clustNames,gapDiff,"ColorScaling","scaledcolumns");
    %     xlabel('No. of Clusters')
    %     h4.title('Gap Value derivative')
    % end
    % 
    % % Plot inconsistency ranges (standard deviation)
    % c=1:1:10;
    % figure('WindowStyle','docked')
    % for c=1:length(c)
    %     f(c) = c;
    %     % Plot inconsistency for lots a range of depths
    %     subplot(5,2,c)
    %     hI(c) = heatmap(linkage_method,distance_method,meani_M(:,:,c))
    %     hI(c).Title = [f(c)]
    % end
    % 
    % % Sample one specific depth for inconsistency value
    % meani_rem = meani_M(:,:,2);
    % 
    % % Heatmaps to visualise results based on cophenet
    % figure('WindowStyle','docked')
    % hC = heatmap(linkage_method,distance_method,c_M)
    % ax = gca;
    % ax.FontSize = 16;
    % 
    % hC.Title = 'Cophenetic correlation';
    % hC.XLabel = 'Linkage Method';
    % hC.YLabel = 'Distance Method';
    % 
    % % Plot heatmap for inconsistency depth (Mean across a range of depths)
    % figure('WindowStyle','docked')
    % hI = heatmap(linkage_method,distance_method,meani_rem)
    % ax = gca;
    % ax.FontSize = 16;
    % hI.Title = 'Inconsistency coefficient';
    % hI.XLabel = 'Linkage Method';
    % hI.YLabel = 'Distance Method';
    
    % Plot violin plots for inconsistency and other clustering metrics
    figure('WindowStyle','docked')
    tl = tiledlayout(6,6,'TileSpacing','tight','Padding','compact');
    set(gcf, 'Renderer', 'painters');
    
    IdxSi = contains(RunName,'single');
    IdxAv = contains(RunName,'average');
    IdxCo = contains(RunName,'complete');
    IdxWa = contains(RunName,'ward');
    
    VioCols = zeros(length(RunName),3);
    VioCols(IdxSi,:) = repmat(IBM_ColorBlind(1,:),sum(IdxSi),1);
    VioCols(IdxAv,:) = repmat(IBM_ColorBlind(2,:),sum(IdxAv),1);
    VioCols(IdxCo,:) = repmat(IBM_ColorBlind(3,:),sum(IdxCo),1);
    VioCols(IdxWa,:) = repmat(IBM_ColorBlind(4,:),sum(IdxWa),1);
    
    % Plot cophenetic correlations
    linkageColors = IBM_ColorBlind(1:length(linkageNames),:);
    
    nEuc = 0;
    nCB = 0;
    nCS = 0;
    for k = 1:length(RunName)
        if contains(RunName{k},'cosine')
            nexttile(1,[2 2])
            nEuc = nEuc+1;
            CC{1}(nEuc) = CopheneticCoeff(k);
            ylabel('Pairwise distance');
        elseif contains(RunName{k},'euclidean')
            nexttile(3,[2 2])
            nCB = nCB+1;
            CC{2}(nCB) = CopheneticCoeff(k);
        elseif contains(RunName{k},'cityblock')
            nexttile(5,[2 2])
            nCS = nCS+1;
            CC{3}(nCS) = CopheneticCoeff(k);
        end
        hold on
    
        if contains(RunName{k},'single'),   col = linkageColors(1,:); Name = 'Single';
        elseif contains(RunName{k},'average'), col = linkageColors(2,:); Name = 'Average';
        elseif contains(RunName{k},'complete'), col = linkageColors(3,:); Name = 'Complete';
        elseif contains(RunName{k},'ward'), col = linkageColors(4,:); Name = 'Ward';
        end
    
        plot(CophenetDists(k,:), PairWiseDists(k,:), '.', 'Color', col);
        MaxDist = max([CophenetDists(k,:), PairWiseDists(k,:)]);
        plot([0 MaxDist],[0 MaxDist],'--k','LineWidth',1);
        box on; set(gca,'layer','top'); axis square;
        xlabel('Cophenetic distance'); 
    end
    
    count = 'a';
    for k = 1:3
        nexttile((k*2)-1,[2 2])
        drawnow
        % Add distance label
        xLimits = xlim;
        yLimits = ylim;
        text(0.9*xLimits(2), ...
             0.02*yLimits(2), ...
             metric_labels{k}, 'FontWeight', 'normal', 'FontSize', 12, ...
             'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
        % Add subplot letter
        text(0.02*xLimits(2), ...
         0.98*yLimits(2), ...
         strcat(count,'.'), 'FontWeight', 'bold', 'FontSize', 12, ...
         'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    
        % Add linkage legend to the each plot
        hold on;
        if k==2
            for i = 1:length(linkageNames)
                lgdplot(i) = plot(nan, nan, '.', 'Color', linkageColors(i,:),'MarkerSize',20, 'DisplayName', string(strcat(linkageNames{i},{' '}, '| CC =',{' '},num2str(round(CC{k}(i),2)))));
            end
        else
            for i = (1:length(linkageNames)-1)
                lgdplot(i) = plot(nan, nan, '.', 'Color', linkageColors(i,:),'MarkerSize',20, 'DisplayName', string(strcat(linkageNames{i},{' '}, '| CC =',{' '},num2str(round(CC{k}(i),2)))));
            end
        end
        %lgdplot(i+1) = plot(nan,nan,'--k','LineWidth',1,'DisplayName','1:1 line');
        lgd = legend(lgdplot, 'Location', 'northeast');
        %title(lgd,'Linkage criteria')
        clear lgd lgdplot
    
        count = char(count+1);
    end
    
    % Inconsistency
    nexttile([2 3])
    meanInconsist = meanInconsist';
    meanInconsist = num2cell(meanInconsist,1);
    vp = daviolinplot(meanInconsist,'violin','half','colors',VioCols,'box',2,'boxcolors','same','outliers',1,'xtlabels',RunName);
    box on
    %title('Inconsistency coefficient')
    ylabel('Inconsistency coefficient')
    hold on
    xline(0.5,'--k','Cosine','LabelHorizontalAlignment','right','LabelVerticalAlignment','bottom','FontWeight','normal','LabelOrientation','horizontal')
    xline(3.5,'--k','Euclidean','LabelHorizontalAlignment','right','LabelVerticalAlignment','bottom','FontWeight','normal','LabelOrientation','horizontal')
    xline(7.5,'--k','Cityblock','LabelHorizontalAlignment','right','LabelVerticalAlignment','bottom','FontWeight','normal','LabelOrientation','horizontal')
    set(gca,'xtick',[])
    xlim([0.5 length(RunName)+0.5]);
    
    xLimits = xlim;
    yLimits = ylim;
    % Add subplot letter
    text(xLimits(1) + 0.02*xLimits(2), ...
     0.98*yLimits(2), ...
     strcat(count,'.'), 'FontWeight', 'bold', 'FontSize', 12, ...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    
    % Silhouette
    count = char(count+1);
    nexttile([2 3])
    siIndex = num2cell(siIndex,1);
    vp = daviolinplot(siIndex,'violin','half','colors',VioCols,'box',2,'boxcolors','same','outliers',1,'xtlabels',RunName);
    box on
    %title('Silhouette coefficient')
    ylabel('Silhouette coefficient')
    hold on
    xline(0.5,'--k','Cosine','LabelHorizontalAlignment','right','LabelVerticalAlignment','bottom','FontWeight','normal','LabelOrientation','horizontal')
    xline(3.5,'--k','Euclidean','LabelHorizontalAlignment','right','LabelVerticalAlignment','bottom','FontWeight','normal','LabelOrientation','horizontal')
    xline(7.5,'--k','Cityblock','LabelHorizontalAlignment','right','LabelVerticalAlignment','bottom','FontWeight','normal','LabelOrientation','horizontal')
    set(gca,'xtick',[])
    xlim([0.5 length(RunName)+0.5]);
    
    xLimits = xlim;
    yLimits = ylim;
    % Add subplot letter
    text(xLimits(1) + 0.02*xLimits(2), ...
     0.98*yLimits(2), ...
     strcat(count,'.'), 'FontWeight', 'bold', 'FontSize', 12, ...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    
    % Legend
    for i = 1:length(linkageNames)
        p(i) = plot(nan,nan,'LineStyle','none','Marker','square','MarkerEdgeColor','k','MarkerFaceColor',IBM_ColorBlind(i,:),'MarkerSize',20,'DisplayName',linkageNames{i});
    end
    lgd2 = legend(p,'Location','northeast');
    
    % Davies Bouldin index
    count = char(count+1);
    nexttile([2 3])
    dbIndex = num2cell(dbIndex,1);
    vp = daviolinplot(dbIndex,'violin','half','colors',VioCols,'box',2,'boxcolors','same','outliers',1,'xtlabels',RunName);
    box on
    %title('Davies Bouldin index')
    ylabel('Davies Bouldin index')
    hold on
    xline(0.5,'--k','Cosine','LabelHorizontalAlignment','right','LabelVerticalAlignment','bottom','FontWeight','normal','LabelOrientation','horizontal')
    xline(3.5,'--k','Euclidean','LabelHorizontalAlignment','right','LabelVerticalAlignment','bottom','FontWeight','normal','LabelOrientation','horizontal')
    xline(7.5,'--k','Cityblock','LabelHorizontalAlignment','right','LabelVerticalAlignment','bottom','FontWeight','normal','LabelOrientation','horizontal')
    set(gca,'xtick',[])
    xlim([0.5 length(RunName)+0.5]);
    
    xLimits = xlim;
    yLimits = ylim;
    % Add subplot letter
    text(xLimits(1) + 0.02*xLimits(2), ...
     0.98*yLimits(2), ...
     strcat(count,'.'), 'FontWeight', 'bold', 'FontSize', 12, ...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    
    % Calinski Harabasz index
    count = char(count+1);
    nexttile([2 3])
    chIndex = num2cell(chIndex,1);
    vp = daviolinplot(chIndex,'violin','half','colors',VioCols,'box',2,'boxcolors','same','outliers',1,'xtlabels',RunName);
    box on
    %title('Calinski Harabasz index')
    ylabel('Calinski Harabasz index')
    hold on
    xline(0.5,'--k','       Cosine','LabelHorizontalAlignment','right','LabelVerticalAlignment','bottom','FontWeight','normal','LabelOrientation','horizontal')
    xline(3.5,'--k','Euclidean','LabelHorizontalAlignment','right','LabelVerticalAlignment','bottom','FontWeight','normal','LabelOrientation','horizontal')
    xline(7.5,'--k','Cityblock','LabelHorizontalAlignment','right','LabelVerticalAlignment','bottom','FontWeight','normal','LabelOrientation','horizontal')
    set(gca,'xtick',[])
    xlim([0.5 length(RunName)+0.5]);
    
    xLimits = xlim;
    yLimits = ylim;
    %ylim([-5 yLimits(2)]); % so Y line label is visible
    % Add subplot letter
    text(xLimits(1) + 0.02*xLimits(2), ...
     0.98*yLimits(2), ...
     strcat(count,'.'), 'FontWeight', 'bold', 'FontSize', 12, ...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

end
