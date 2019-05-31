function sol=CA_sol_double(cells,conn,lambda,time,stop_stim,...
	t_r,seed,n_conn,overlap)
% Cellular automaton where cells with 'n_conn' or more neighbors
% can only fire when stimulated or two neighbors are firing at the
% same time. The CA is run from time t=1 to t='time' and the
% results are returned in 'sol', where the dimensions of 'sol' are
% [length(cells) time].
% 'cells' gives the initial conditions in an array.
% 'lambda' is the mean stimulation rate for the Poisson process
% that governs random external stimulation.
% 'stop_stim' is the last time step for stimulation.
% 't_r' is the refractory period for all cells.
% 'seed' is the random seed for the Poisson process.
% 'overlap' gives "firing" length. If overlap=3,
% a cell is "firing" when it is exctied (for one time
% step), and then the next 3 time steps.

n_cells = length(cells);
sol=zeros(n_cells,time);
sol(:,1)=cells;

rand('seed',seed);
stim = zeros(n_cells,time);
p = lambda*exp(-lambda);
stim(:,1:stop_stim) = rand(n_cells,stop_stim)<p;

for k=2:time
	for i=1:n_cells
		if sol(i,k-1)==1
			sol(i,k)=-t_r;
		elseif sol(i,k-1)<0
			sol(i,k) = sol(i,k-1)+1;
		else % sol(i,k)==0
			if next_to_excited(i,sol(:,k-1),conn,t_r,n_conn,overlap)...
					|| stim(i,k-1)
				sol(i,k) = 1;
			end
		end
	end
end

% subfunctions-----------------------------
function yes=next_to_excited(index,sol_prev,conn,...
	t_r,n_conn,overlap)

yes=0;
c=conn{index};
num_exc = 0;
num_overlap = 0;
for j=1:length(c)
	if sol_prev(c(j))==1 
		num_exc = num_exc + 1;
	elseif sol_prev(c(j))<-t_r+overlap
		num_overlap = num_overlap + 1;
	end
end
if length(c)<n_conn && num_exc>0
	yes = 1;
elseif num_exc>=1 && num_exc+num_overlap>1
	yes = 1;
end


