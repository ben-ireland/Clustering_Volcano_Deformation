function NumClusts = recolorDendrogramClasses(H, T, labels, Cmap, Lab)
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
% recolorDendrogramClasses
% Recolors a dendrogram based on user-defined cluster assignments.
%
% Inputs:
%   H      - handles returned by dendrogram/dendrogram2
%   T      - cluster assignments (vector of length N), 1..K
%   labels - leaf labels used in original dendrogram call
%   Cmap   - K×3 RGB colormap. Row i = color of cluster i.
%
% Behavior:
%   - All dendrogram branches are recolored based on T.
%   - Branches that span >1 cluster remain black.
%   - Supports left or top orientation.
%   - Works with numeric or string labels.

    % --- Extract leaf order from axis tick labels ---
    ax = ancestor(H(1), 'axes');
    xlabels = get(ax, 'XTickLabel');
    ylabels = get(ax, 'YTickLabel');

    % Determine dendrogram orientation:
    orientationLeft = length(ylabels) - length(H) ==1;  

    % Handle both orientations
    if orientationLeft
        leafLabels = ylabels;    % leaves are on y-axis
    else
        leafLabels = xlabels;    % default dendrogram
    end

    % Convert labels to numeric indicies
    for i = 1:length(leafLabels)
        obsIndex(i) = str2double(leafLabels(i,:));
    end

    % --- Recolor each dendrogram line ---
    numCol = 0;
    numBlk = 0;
    for i = 1:length(H)
        xd = get(H(i), 'XData');
        yd = get(H(i), 'YData');

        % Get leaf positions on the axis
        if orientationLeft
            yvals = yd(~isnan(yd));
        else
            xvals = xd(~isnan(xd));
        end

        leafPositions = unique(round(yvals));

        % Convert to actual observation indices
        leaves = obsIndex(leafPositions);

        % Get cluster IDs of these leaves
        leafClusters = unique(T(leaves));

        if numel(leafClusters) == 1
            numCol = numCol + 1;
            % Pure branch—assign color for this cluster:
            cID = leafClusters(1);
            cIDcount(numCol) = leafClusters(1);
            LineColCount(numCol) = i;
            %set(H(i), 'Color', Cmap(cID, :));
        else
            numBlk = numBlk + 1;
            LineBlackCount(numBlk) = i;
            % Mixed branch—set black
            %set(H(i), 'Color', [0 0 0]);
        end
    end
    
    % Correct colourmap index sequentially when values have been removed
    % because they are black
    cIDNums = unique(cIDcount);
    cIDNumsCorrect = 1:length(cIDNums);
    cIDcount2 = zeros(size(cIDcount));
    NumClusts = length(cIDNums);
    for i = 1:length(cIDNums)
        % Re-number cIDs to be sequential
        idxThiscID = cIDcount == cIDNums(i);
        cIDcount2(idxThiscID) = cIDNumsCorrect(i);
    end
    cIDcount = cIDcount2;
    
    for i = 1:length(LineColCount)
        set(H(LineColCount(i)), 'Color', Cmap(cIDcount(i), :));
    end

    set(H(LineBlackCount), 'Color', [0 0 0]);

    % Add in Labels
    if Lab ==1
        if orientationLeft
            set(gca,'YTickLabel',labels(obsIndex));
        end
    end
end
