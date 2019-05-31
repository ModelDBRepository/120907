function cycles = find_cycles(conn,node)
% Give histogram of cycle lengths for each neighbor of 'node'.
% Cycles is a cell array, where entry i is a 2 column matrix for
% neighbor i. Each row of the two column matrix gives
% [cycle_length #cycles_with_that_length]. Note that each cycle is
% accounted for twice, since it must go through two neighbors of 'node'.

lc = length(conn);

% Make a spanning tree centered around node
% where each vertex in the tree is placed so that
% its shortest distance in the graph is its distance to node.
% Keep previous vertex
% Keep address for each vertex based on branch points.
% Keep distance from each vertex to closest branch point.
% Flag for branch when just AFTER branch node to tell
%    which direction to go from tree

E_tree = []; % edges in tree
E_cycle = []; % edges that make a cycle
% 1st_address1 1st_address2 cell1 cell2

rank_tree = repmat(struct('prev',[],'address',[],...
	'dist',[],'branch',[]),lc,1);
% 'prev' = previous vertex
% 'address' = array of vertices just beyond branch points, starting
% with immediate neighbor of 'node'
% 'dist' = distance from vertex to closest branch point
% 'branch' = is this vertex just beyond a branch point?

S = []; % searched
R = node; % reached
rank_tree(node) = struct('prev',0,'address',[],...
	'dist',0,'branch',0);

while ~isempty(R)
  %fprintf(1,'current cell: %d\n',R(1));
  nc = length(conn{R(1)});
  % decide if we are at branching point,
  % make note to add address if we are
  if nc-length(intersect([R S],conn{R(1)}))>=2
    % if the number of neighbors not yet reached is at least 2
    branch_node = 1;
  else
    branch_node = 0;
  end
  
  for k=1:length(conn{R(1)})
    if isempty(find([R S]==conn{R(1)}(k),1))
      %fprintf(1,'not reached yet: %d\n',conn{R(1)}(k));
      R = [R conn{R(1)}(k)];
      rank_tree(conn{R(1)}(k)).prev = R(1);
      
      % is R(1) part of the address, 
      % i.e. branch of previous branch point?
      if rank_tree(R(1)).branch
        rank_tree(conn{R(1)}(k)).dist = 1;
        rank_tree(conn{R(1)}(k)).address = ...
            [rank_tree(R(1)).address R(1)];
      else
        rank_tree(conn{R(1)}(k)).dist = ...
            rank_tree(R(1)).dist + 1;
        rank_tree(conn{R(1)}(k)).address = ...
            rank_tree(R(1)).address;
      end
      if branch_node
        rank_tree(conn{R(1)}(k)).branch = 1;
      else
        rank_tree(conn{R(1)}(k)).branch = 0;
      end
    elseif conn{R(1)}(k) ~= rank_tree(R(1)).prev
      %fprintf(1,'already in tree: %d\n',conn{R(1)}(k));
      if R(1)<=conn{R(1)}(k) && (isempty(E_cycle) || ...
                                 ~any(E_cycle(:,3)==R(1) & ...
                                      E_cycle(:,4)==conn{R(1)}(k)))
        E_cycle = [E_cycle; branch(R(1),rank_tree)...
                   branch(conn{R(1)}(k),rank_tree)...
                   R(1) conn{R(1)}(k)];
      elseif isempty(E_cycle) || ...
            ~any(E_cycle(:,4)==R(1) &...
                 E_cycle(:,3)==conn{R(1)}(k))
        E_cycle = [E_cycle; ...
                   branch(conn{R(1)}(k),rank_tree)...
                   branch(R(1),rank_tree)...
                   conn{R(1)}(k) R(1)];
      end
    end
  end
  S = [S R(1)];
  R = R(2:length(R));
end

% Cycles going through node are made from combinations
% of edges in E_cycle. Just compute cycles with 
% one extra edge
n_branches = length(conn{node});
sE = size(E_cycle);
cycles = cell(n_branches,1);

for i=1:sE(1)
  if E_cycle(i,1)==E_cycle(i,2)
    % cycle is on same branch, doesn't go through node
    continue;
  end
  
  % find length of cycle through E_cycle(i,:)
  for j=3:4
    dist(j) = rank_tree(E_cycle(i,j)).dist;
    for k=1:length(rank_tree(E_cycle(i,j)).address);
      stop = rank_tree(E_cycle(i,j)).address(k);
      dist(j) = dist(j)+rank_tree(stop).dist;
    end
  end
  cycle_length = sum(dist)+1;
  
  % add cycle to histogram
  for j=1:2
    branch_ind = find(conn{node}==E_cycle(i,j),1);
    sc = size(cycles{branch_ind});
    if sc(1)>0
      cycle_ind = find(cycles{branch_ind}(:,1)...
                       ==cycle_length,1);
    else
      cycle_ind = [];
    end
    if isempty(cycle_ind)
      cycles{branch_ind} = [cycles{branch_ind};...
                          cycle_length 1];
    else
      cycles{branch_ind}(cycle_ind,2) = ...
          cycles{branch_ind}(cycle_ind,2)+1;
    end
  end
end


function b=branch(cell,rank_tree)
% give branch of cell based on address

if isempty(rank_tree(cell).address)
  b = cell;
else
  b=rank_tree(cell).address(1);
end
