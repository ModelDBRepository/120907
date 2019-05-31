function cluster = conn_analysis(conn)
% Analyze graph given by cell array conn. Print out a list of
% network properties and return a cell array of
% clusters, where cluster{i} is an array of all cells in cluster
% i. Clusters are simply listed in order that they are found
% according to the algorithm. Larger clusters tend to be found first.

sc = length(conn); % number of cells
fprintf(1,'Total number of cells: %d\n',sc);

% count number of connections each cell has
nempty = 0;
n = zeros(1,4);
n_invalid = 0;
for i=1:sc
	if isempty(conn{i})
		nempty = nempty + 1;
	elseif length(conn{i}) <= 4
		lg = length(conn{i});
		n(lg) = n(lg) + 1;
	else
		n_invalid = n_invalid + 1;
	end
end
fprintf(1,'Number of cells connected to\n')
fprintf(1,'0 other cells: %d\n',nempty)
fprintf(1,'1 other cell: %d\n',n(1))
fprintf(1,'2 other cells: %d\n',n(2))
fprintf(1,'3 other cells: %d\n',n(3))
fprintf(1,'4 other cells: %d\n',n(4))
fprintf(1,'>4 other cells: %d\n',n_invalid)

% Check that connections are reciprocal
for i=1:sc
	for j=1:length(conn{i})
		if isempty(find(conn{conn{i}(j)}==i,1))
			fprintf(1,'Cell %d is not reciprocal\n',i);
		end
	end
end

% find clusters
cells = 1:sc; % list of cells not yet dealt with
j = 1; % current cluster
while ~isempty(cells)
	cells_left = length(cells);
	cell = cells(1); % take first cell not dealt with
	cells = cells(2:length(cells)); %remove it from list
	cluster{j} = cell; % make a new cluster with it
	ends = conn{cell}; % get its connections
  % keep following connections until every cell in cluster is added
  % to cluster{j}
	while ~isempty(ends)
		cluster{j} = [cluster{j} ends]; %add cells to cluster
		% remove ends from cells one at a time
		for k = ends 
			ki = find(cells == k);
			if ki > 1
				front = 1:(ki-1);
			else front = []; end
			if ki < length(cells)
				back = (ki+1):length(cells);
			else back = []; end
			cells = [cells(front) cells(back)];
			cells_left = length(cells);
		end
		old_ends = ends;
		ends = [];
		for k = old_ends % for all old ends
			% test each cell they're connected to
			for test = conn{k}
				% to see if they're in the cluster or ends
				inclus = find(test == cluster{j},1);
				inends = find(test == ends,1);
				for ci = 1:j-1
					if find(test == cluster{ci})
						error('end in other cluster')
					end
				end
				if isempty(inclus) && isempty(inends)
					ends = [ends test];
				end
			end
		end % k = old_ends
	end % ~isempty(ends)
	% now that all cells in the cluster are added
	% move on to new cluster
	cells_left = length(cells);
	j = j + 1;
end % ~isempty(cells)

fprintf(1,'The number of clusters is: %d\n',length(cluster));

% count the number of cells in each cluster
cluster_count = [length(cluster{1}); 1];
for i = 2:length(cluster)
	j = find(length(cluster{i})==cluster_count(1,:));
	if isempty(j)
		new_col = [length(cluster{i}); 1];
		cluster_count = [cluster_count new_col];
	else
		cluster_count(2,j) = cluster_count(2,j)+1;
	end
end
scc = size(cluster_count);
for i=1:scc(2)
	fprintf(1,'There are %d clusters with %d cells.\n',...
		cluster_count(2,i),cluster_count(1,i));
end
