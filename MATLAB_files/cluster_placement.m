function pm = cluster_placement(clusters,cells)
% Take the cell array clusters which contains
% which cells are in which clusters, and order the
% cells from smallest cluster to largest cluster.
% Return pm where pm(i) gives where cells(i) is in order.
% This assumes that cells is in ascending order.

lc = length(clusters);
cluster_lengths = zeros(1,lc);
for i = 1:lc
	cluster_lengths(i) = length(clusters{i});
end

[srtd_clust srtd_ind] = sort(cluster_lengths);

cell_list = [];
for i = 1:lc
	for j = clusters{srtd_ind(i)}
		if ~isempty(find(cells==j,1))
			cell_list = [cell_list j];
		end
	end
end

[c pm] = sort(cell_list);