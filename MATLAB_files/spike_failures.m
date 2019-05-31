function [X Y] = spike_failures(t,state,cells,conn)
% Give the times that cells fail to fire even though
% one of their neighbors fires.
% X for failures, Y for spikes
% X{1} for cells with 0 connections
% X{2} for cells with 1 connections
% X{3} for cells with 2 connections
% X{4} for cells with 3 connections
% X{5} for cells with 4 connections
% Row i for X gives [cells(i) 'time of spike failure']
% Row i for Y gives [cells(i) 'time of spike start']
% The time of spike start/failure is when the voltage crosses 5 mV
% from below. Is the voltage reaches over 100 mV, then it is
% counted as a spike, otherwise it is a failure.

X = cell(5,1); 
Y = cell(5,1);
for i=1:length(cells)
	spike_start = find(state(:,cells(i))>=5,1)
	while ~isempty(spike_start)
        spike_length = ...  
					find(state(spike_start+1:length(t),cells(i))...
					<5,1);
				n_conn = length(conn{cells(i)});
		if max(state(spike_start:...    
				spike_start+spike_length,cells(i)))...
				<100
			X{n_conn+1} = ...
				[X{n_conn+1}; cells(i) t(spike_start)];
		else
			Y{n_conn+1} = ...
				[Y{n_conn+1}; cells(i) t(spike_start)];
		end
		spike_start = ...
			find(state(spike_start+spike_length+1:...
			length(t),cells(i))>=5,1)...
			+ spike_start+spike_length
	end
end