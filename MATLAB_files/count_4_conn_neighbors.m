function [distribution hh] = count_4_conn_neighbors(...
	cells,conn)
% This routine is used to determine what is a "typical"
% gap-junctional load for the small network simulations. In order
% to do this, for each 4-connected cell in the large cluster (the
% array 'cells') we took turns choosing which neighbor would be
% "upstream" and then counted how many connections each
% "downstream" neighbor had. The total "downstream" connections is
% an approximation of the gap-junctional load.
% 'cells' lists the cells that we want to take into account.
% 'conn' gives the cell connections.
% 'distribution' is an array where distribution(i) is the number of
% instances found where a 4-connected cells has i-1 down stream
% second tier neighbors.
% 'hh' is an array of extra heavily connected cells, number of
% instances >= 8.

hh = [];

distribution = zeros(1,10);
for i=1:length(cells)
	if length(conn{cells(i)})==4
		for j=1:4
			% for each combo without neighbor j
			sum = 0;
			for k=1:4
				if k~=j
					sum = sum+length(conn{conn{cells(i)}(k)})-1;
				end
			end
			distribution(sum+1) = distribution(sum+1)+1;
			if sum>=8
				hh = [hh i];
			end
		end
	end
end