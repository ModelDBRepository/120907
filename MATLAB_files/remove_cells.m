function sconn = remove_cells(conn,cells)
% Remove all connections of each cell listed in 'cells' and return
% the resulting graph in 'sconn'. Original connections are
% specified in 'conn'.

sconn = conn;
rem_list = [];
for i=1:length(sconn)
	if ~isempty(find(cells==i))
		rem_list = [rem_list i];
	end
end

for i=rem_list
	% remove connections to cell
	for j = 1:length(sconn{i})
		oc = sconn{sconn{i}(j)};
		k = find(oc == i);
		sconn{sconn{i}(j)} = [oc(1:k-1) oc(k+1:length(oc))];
	end
	% remove connections
	sconn{i} = [];	
end
