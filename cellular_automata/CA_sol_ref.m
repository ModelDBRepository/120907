function sol=CA_sol_ref(cells,conn,lambda,time,t_r,seed,...
	stop_stim)
% Cellular automata where 4-connected cells have a longer
% refractory period than cells with fewer connections. The CA is
% run from time t=1 to t='time' and the results are returned in
% 'sol', where dimensions of 'sol' are = [length(cells) time].
% 'cells' gives the initial conditions in an array.
% 'lambda' is the mean stimulation rate of the Poisson process that
% governs random external stimulation.
% Cells with 4 connections have a refractory period
% of t_r(2) instead of t_r(1).
% 'seed' gives random seed to base simulation on. 
% 'stop_stim' gives time to stop the stimulation.
% If 'stop_stim'=='time'. there is stimulation whole time.
% Assume that 'stop_stim' <= 'time'.

rand('seed',seed);
n_cells = length(cells);

sol=zeros(n_cells,time);
sol(:,1)=cells;
stim = zeros(n_cells,time);
p = lambda*exp(-lambda); % probability of stim
stim(:,1:stop_stim) = rand(n_cells,stop_stim)<p;

for k=2:time
	for i=1:n_cells
		if sol(i,k-1)==1
			if length(conn{i})==4
				sol(i,k) = -t_r(2);
			else
				sol(i,k) = -t_r(1);
			end
		elseif sol(i,k-1)<0
			sol(i,k) = sol(i,k-1)+1;
		elseif next_to_excited(i,sol(:,k-1),conn) ...
				|| stim(i,k-1)
			sol(i,k) = 1;
		end
	end
end

% subfunctions-----------------------------
function yes=next_to_excited(index,sol_prev,conn)

yes=0;
c=conn{index};
for j=1:length(c)
	if sol_prev(c(j))==1
		yes=1;
		break;
	end
end

