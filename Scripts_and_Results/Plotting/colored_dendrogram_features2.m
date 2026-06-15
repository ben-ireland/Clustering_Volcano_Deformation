function colored_dendrogram_features2(Z, X, L, N, C, T, TS)
%% Ben Ireland, June 2026, Cluster multi-variate datasets
% School of Earth Sciences, University of Bristol
%
% Methods in these scripts are described in:
% "Clustering global volcano deformation
% datasets: insights and limitations for analogue signals"
% Submitted to JVGR
%
% Copyright: Ben Ireland, 2026
% V1.0
%
% COLORED_DENDROGRAM_FEATURES Draw dendrograms with colored horizontal & vertical lines
% reflecting the mean of a feature across subclusters.
%
%   colored_dendrogram_features(Z, X)
%
% Inputs:
%   Z - linkage matrix (as from linkage())
%   X - (n x p) matrix, rows = observations, columns = features
%   L - Feature labels (p x 1 cell array of strings)
%   N - Desired number of clusters e.g. 10
%   C - Discrete colourmap for clusters (N x 3 matrix)
%   T - Overall title for plot (string)
%   TS - Plot overall title (1) or not (0)
%
% For each feature, plots a dendrogram where:
%   • Each vertical line is colored by the mean feature value of its subcluster
%   • Each horizontal link is colored by the mean feature value of the merged cluster
%
%
% -------------------------------------------------------------------------

    [n, p] = size(X);
    m = size(Z,1);
    totalClusters = n + m;

    % --- Compute average feature values per cluster ---
    clusterVals = zeros(totalClusters, p);
    clusterVals(1:n,:) = X;
    clusterSizes = ones(totalClusters,1);

    for i = 1:m
        c1 = Z(i,1);
        c2 = Z(i,2);
        idx = n + i;
        s1 = clusterSizes(c1);
        s2 = clusterSizes(c2);
        clusterVals(idx,:) = (s1*clusterVals(c1,:) + s2*clusterVals(c2,:)) / (s1+s2);
        clusterSizes(idx) = s1 + s2;
    end

    % --- Get dendrogram structure (without plotting lines) ---
    [~,~,perm] = dendrogram(Z, 0);
    close(gcf); % remove temporary figure

    % --- Compute x-positions for each cluster in leaf order ---
    xPos = compute_x_positions(Z, perm);

    % --- Plot one dendrogram per feature ---
    figure('WindowStyle','docked')
    t = tiledlayout(ceil((p+1)/2),2, 'TileSpacing','tight', 'Padding','compact');
    
    if TS ==1
        title(t,T,'Interpreter','none')
    end
    count = 'a';

    nexttile
    cutoff_for_colours = median([Z(end-N+1,3) Z(end-N+2, 3)]);
    [H, ~, P, ~,theGroups] = dendrogram2(Z, 0, 'Labels',[], 'orientation','left','ColorThreshold',cutoff_for_colours);
    hold on
    set(H,LineWidth=2)

    T = cluster(Z,"maxclust",N);                 % your cluster() output
    Cmap = C(2:end,:);
    labels = L;
    NumMultClusts = recolorDendrogramClasses(H, T, labels, Cmap, 0);

    % Add dummy legend to plot 1
    for j = 1:NumMultClusts
        hold on
        dendLgd(j) = plot(nan,nan,'-','Color',Cmap(j,:),'LineWidth',2,'DisplayName',['Cluster ', num2str(j)]);
    end

    if NumMultClusts~=N
        dendLgd(end+1) = plot(nan,nan,'-','Color','k','LineWidth',2,'DisplayName','Outlier');
    end

    dendLgd(end+1) = xline(cutoff_for_colours,':','LabelHorizontalAlignment','left','LabelVerticalAlignment','top','LabelOrientation','horizontal','DisplayName','Cut-off','LineWidth',2,'Color',[0.5 0.5 0.5]);
    lgdDend = legend(dendLgd,'Location','southwest');
    
    
    hold on
    title('Clusters')
    box on; set(gca,'Layer','top');
    if TS==1
        xlabel('Linkage Height')
    end

    yyaxis right
    set(gca,'YTick',[]);

    yyaxis left % Format axes same as the rest of the charts
    Lim = ceil(size(T,1)/10);
    LimMax = Lim*10;
    ylim([0 LimMax])
    set(gca,'YTick',0:50:LimMax);
    set(gca,'YTickLabel',string(0:50:LimMax));
    if TS ==0
        set(gca,'YTick',[]);
    else
        ylabel('Observation order')
    end

    ax = gca;
    ax.YAxis(1).Color = [0 0 0];
    ax.YAxis(2).Color = [0 0 0];

%     Plt = xline(cutoff_for_colours,':','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','DisplayName','Cut-off value','Color',[0.5 0.5 0.5],'LineWidth',2);
%     % Dummy legend
%     for k = 1:size(colors,1)
%         if k ~=size(colors,1)
%             Plt2(k) = plot(nan,nan,'-',LineWidth=2,Color=colors(k,:),DisplayName=['Cluster ',num2str(k)]);
%         else
%             Plt2(k) = plot(nan,nan,'-',LineWidth=2,Color=colors(k,:),DisplayName='Outlier');
%         end
%     end
% 
%     lgd = legend([Plt2 Plt],'Location','southwest');

    xLimits = xlim;
    yLimits = ylim;

    text(0.99*xLimits(2), ...
     0.99*yLimits(2), ...
     strcat(count,'.'), 'FontWeight', 'bold', 'FontSize', 12, ...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

    % First dendrogram (normal colours for clusters)
    for f = 1:p
        count = char(count+1);

        nexttile; hold on;
        title(L{f});
        if TS==1
            xlabel('Observation order');
            ylabel('Linkage height');
        else
            set(gca,'XTick',[]);
            if f==p || f==p-1
                ylabel('Linkage height')
            end
        end
        set(gca,'YDir','normal');

        vals = clusterVals(:,f);
        maxval = max(abs(vals));

        if min(X(:,f))<0
            colormap(gca, flipud(cbrewer2('RdYlBu',256)));
            cmap = flipud(cbrewer2('RdYlBu',256));
            clim = [-maxval maxval];
        else
            colormap(gca, cbrewer2('YlOrRd',256));
            cmap = cbrewer2('YlOrRd',256);
            clim = [min(vals) max(vals)];
        end
        

        % --- Draw vertical & horizontal lines for each merge ---
        for i = 1:m
            c1 = Z(i,1);
            c2 = Z(i,2);
            h = Z(i,3);
            xp1 = xPos(c1);
            xp2 = xPos(c2);
            xc  = mean([xp1 xp2]);
            hp1 = height_of_cluster(c1, Z, n);
            hp2 = height_of_cluster(c2, Z, n);

            % Values
            v_parent = clusterVals(n+i,f);
            v_c1 = clusterVals(c1,f);
            v_c2 = clusterVals(c2,f);

            % Colors
            c_parent = val2color(v_parent, clim, cmap);
            c1_col = val2color(v_c1, clim, cmap);
            c2_col = val2color(v_c2, clim, cmap);

            % --- Vertical lines ---
            plot([xp1 xp1], [hp1 h], 'Color', c1_col, 'LineWidth', 2);
            plot([xp2 xp2], [hp2 h], 'Color', c2_col, 'LineWidth', 2);

            % --- Horizontal line connecting them ---
            plot([xp1 xp2], [h h], 'Color', c_parent, 'LineWidth', 2);
        end

        caxis(clim);
        cb = colorbar;
        ylabel(cb, 'Feature value');
        Plt = yline(cutoff_for_colours,':','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','DisplayName','Cut-off value','Color',[0.5 0.5 0.5],'LineWidth',2);
        
        view([-90 90])

        box on
        set(gca,'Layer','top')

        xlim(yLimits)
        ylim(xLimits);
        text(0.99*yLimits(2), ...
         0.99*xLimits(2), ...
         strcat(count,'.'), 'FontWeight', 'bold', 'FontSize', 12, ...
         'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
    end
end


% -------------------------------------------------------------------------
% Helper: Map value to RGB color from a colormap
function c = val2color(v, clim, cmap)
    nC = size(cmap,1);
    t = (v - clim(1)) / (clim(2)-clim(1));
    t = max(0, min(1, t));
    idx = max(1, round(t*(nC-1))) + 1;
    c = cmap(idx,:);
end

% -------------------------------------------------------------------------
% Helper: Compute x-positions for each cluster in dendrogram
function xPos = compute_x_positions(Z, perm)
    n = length(perm);
    m = size(Z,1);
    total = n + m;
    xPos = zeros(total,1);
    % Leaves get fixed positions
    xPos(perm) = 1:n;
    % Internal nodes: mean of children
    for i = 1:m
        c1 = Z(i,1);
        c2 = Z(i,2);
        xPos(n+i) = mean([xPos(c1) xPos(c2)]);
    end
end

% -------------------------------------------------------------------------
% Helper: Get the height (linkage distance) of a cluster
function h = height_of_cluster(idx, Z, n)
    if idx <= n
        h = 0;
    else
        h = Z(idx - n, 3);
    end
end
