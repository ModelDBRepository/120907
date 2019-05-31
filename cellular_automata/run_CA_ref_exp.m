function run_CA_ref_exp
% Run CA_sol_ref with different refractory periods. Cells with less
% than 4 connections have t_r = 10, which 4-connected cells have
% t_r = 12 + 0:12. Run 10 simulations for each different t_r for
% 4-connected cells. Save output in '.mat' files.

t_r = 10;
dtr = 0:12;
seed = 1:10;
time = 400;
lambda = .0005;
% lambda ~ 2/s/cell, 1 time step ~ 0.25 ms

cells = zeros(32*96,1);
conn = load('conn.out');
conn = cconn2mconn(conn);
%conn = connections(32,96,10,1.6);
sol = cell(length(seed),1);

for i=1:length(dtr)
	for j=1:length(seed)
		sol{j} = CA_sol_ref(cells,conn,lambda,time,...
			[t_r t_r+dtr(i)],seed(j),200);
		file = sprintf('CA_ref_exp_%d.mat',dtr(i));
		dt_r = dtr(i);
		save(file,'t_r','dt_r','seed','sol','lambda')
	end
end


	