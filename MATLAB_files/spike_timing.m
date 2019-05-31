function X = spike_timing(t,state,placement,stim,cells,p)
% Calculate the spike times of every neuron in the network, and
% plot a rastergram.
% 't' is a column vector listing every time step
% 'state' is a matrix, where row i lists the voltage of every cell
% at time t(i)
% 'placement' is an array specifying the order the cells should be
% plotted in - placement(i) gives position of cells(i)
% 'stim' is specifies stimulation to cells where each row gives [cell
% stim_time]. Note that since this data is usually taken from C
% code, the first cell is 0 instead of 1.
% 'cells' is an array specifying which cells to plot: for instance,
% the data for cell 3 is in the third column of 'state'
% 'p' is a boolean indicating whether to plot or not
% 'X' is a 2 column matrix where row i gives ['placement of cell i'
% 'time of spike']

ss = size(state);
lt = length(t);
j=0;
X = [];
for i=1:length(cells) % each neuron
	ind = find(state(:,cells(i))>=50,1);
	while length(ind)>0;
		up = ind(1);
		ind = find(state((up+1):length(t),cells(i))<=50,1);
		if length(ind)>0
			down = ind(1);
		else
			down = 0;
		end
		ti = up + round(down/2);
		X = [X; placement(i) t(ti)];
		j = j+1;
		ind = find(state(up+down+1:ss(1),cells(i))>=50,1)+up+down;
	end
end

if p==1
	plot(X(:,2),X(:,1),'.r');
	hold on
	plot_stim = [];
	ss = size(stim);
	for i = 1:ss(1)
		ind = find(cells==(stim(i,1)+1));
		if ~isempty(ind) && stim(i,2)>=t(1) && stim(i,2)<=t(lt)
			plot_stim = [plot_stim;ind stim(i,2)];
		end
	end
	sp = size(plot_stim);
	if sp(1)>0
		plot(plot_stim(:,2),placement(plot_stim(:,1)),'.b');
	end
	hold off
end

 h = findobj(gca,'Type','line')
 set(h,'MarkerSize',3)