function [active inactive] = heavily_connected(...
	conn,cells,state,soma)
% Give list of 4-connected cells.
% 'conn' is a cell array specifying the connections
% 'cells' is the array of cells to determine whether they're
% heavily connected.
% 'state' is matrix with voltage of 'cells'; the voltage of cell i
% is in column i of 'state'.
% 'soma' indicates whether these are somatic voltages (cutoff 1 mV
% for 'active', or axonal voltage (cutof 50 mV for 'active')
% 'active' is an array of 4-connected cells that fired at least
% once.
% 'inactive' is an array of 4-connected cells that never fired.

active = [];
inactive = [];
for i=1:length(cells)
	if length(conn{cells(i)}) >= 4
		if soma
			if ~isempty(find(state(:,i)>=1,1))
				active = [active i];
			else
				inactive = [inactive i];
			end
		else
			if ~isempty(find(state(:,i)>=50,1))
				active = [active i];
			else
				inactive = [inactive i];
			end
		end
	end
end