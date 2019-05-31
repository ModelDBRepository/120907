function run_CA_double_exp
% Run CA_sol_double with different with 10 different seeds, as well
% as different values of 'overlap' and 'n_conn'. Save output in a .mat file

t_r = 10;
overlap = 0:3;
n_conn = 3:4;
seed = 1:10;
time = 400;
lambda = .0005;
% lambda ~ 2/s/cell, 1 time step ~ 0.25 ms

cells = zeros(32*96,1);
conn = load('conn.out');
conn = cconn2mconn(conn);
%conn = connections(32,96,10,1.6);
sol = cell(length(overlap),length(n_conn),length(seed));

for i=1:length(overlap)
	for j=1:length(n_conn)
		for k=1:length(seed)
			sol{i,j,k} = CA_sol_double(cells,conn,lambda,time,...
				200,12,seed(k),n_conn(j),overlap(i));
			save('CA_double_exp.mat','t_r','seed','sol',...
				'overlap','n_conn','lambda')
		end
	end
end

