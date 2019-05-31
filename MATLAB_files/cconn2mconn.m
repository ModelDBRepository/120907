function mconn = cconn2mconn(cconn)
% Convert conn.out matrix from c output to cell array readable to
% matlab routines. C output is given in an Ncell * 4 matrix. Row i
% lists the cells connected to cell i, with a maximum of 4
% connections. If cell i has less than 4 connections, extra spots
% in the row are filled with -1.

sc = size(cconn);

for i=1:sc(1)
	j = find(cconn(i,:)==-1,1);
	if length(j) ==0
		j = 5;
	end
	mconn{i} = cconn(i,1:(j-1)) + 1;
end