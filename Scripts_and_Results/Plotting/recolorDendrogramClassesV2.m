function NumClusts = recolorDendrogramClassesV2(H,Z,T,Cmap,Lab,labels)
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
% recolorDendrogramClassesV2
%
% Recolour dendrogram branches using exact linkage-tree membership.
%
% Inputs
% -------
% H       : handles returned by dendrogram
% Z       : linkage matrix
% T       : cluster assignments (N x 1)
% Cmap    : K x 3 colour map
% Lab     : relabel leaves flag
% labels  : leaf labels
%
% Output
% -------
% NumClusts : number of represented clusters

NumClusts = 0;

if isempty(H)
    return
end

T = T(:);

N = size(Z,1)+1;

if length(T) ~= N
    error('Length(T) must equal number of observations in Z.');
end

ax = ancestor(H(1),'axes');

%% --------------------------------------------------------------------
% Build exact descendant list for every node
%% --------------------------------------------------------------------

nodeLeaves = cell(2*N-1,1);

for i = 1:N
    nodeLeaves{i} = i;
end

pureNode = false(N-1,1);
nodeCluster = nan(N-1,1);

for r = 1:(N-1)

    left  = Z(r,1);
    right = Z(r,2);

    leaves = [ ...
        nodeLeaves{left}, ...
        nodeLeaves{right} ];

    nodeLeaves{N+r} = leaves;

    clusts = unique(T(leaves));

    if numel(clusts)==1

        pureNode(r) = true;
        nodeCluster(r) = clusts;

    end
end

%% --------------------------------------------------------------------
% Determine dendrogram orientation
%% --------------------------------------------------------------------

xrange = diff(get(ax,'XLim'));
yrange = diff(get(ax,'YLim'));

orientationLeft = xrange < yrange;

%% --------------------------------------------------------------------
% Determine merge height for each line handle
%% --------------------------------------------------------------------

lineHeight = nan(length(H),1);

for i = 1:length(H)

    xd = get(H(i),'XData');
    yd = get(H(i),'YData');

    if orientationLeft

        lineHeight(i) = max(xd(~isnan(xd)));

    else

        lineHeight(i) = max(yd(~isnan(yd)));

    end
end

%% --------------------------------------------------------------------
% Match line handles to linkage rows
%% --------------------------------------------------------------------

mergeHeight = Z(:,3);

usedLines = false(size(H));

for r = 1:(N-1)

    h = mergeHeight(r);

    idx = find( ...
        abs(lineHeight - h) < 1e-10 & ...
        ~usedLines );

    if isempty(idx)
        continue
    end

    idx = idx(1);

    usedLines(idx) = true;

    if pureNode(r)

        cID = nodeCluster(r);

        if cID <= size(Cmap,1)

            set(H(idx), ...
                'Color', Cmap(cID,:));

        else

            set(H(idx), ...
                'Color', [0 0 0]);

        end

    else

        set(H(idx), ...
            'Color', [0 0 0]);

    end
end

%% --------------------------------------------------------------------
% Count represented clusters
%% --------------------------------------------------------------------

NumClusts = numel(unique(nodeCluster(~isnan(nodeCluster))));

%% --------------------------------------------------------------------
% Optional relabelling
%% --------------------------------------------------------------------

if nargin >= 6 && Lab

    try

        if orientationLeft

            set(ax,'YTickLabel',labels);

        else

            set(ax,'XTickLabel',labels);

        end

    catch

        warning('Unable to relabel dendrogram.');

    end

end

end