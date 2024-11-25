function [dA_flipped_connected_sorted, swc_data_updated] = process_floating_branches(dA, swc_data, true_root_id, root_node_ids, parent_node_ids)
    % PROCESS_FLOATING_BRANCHES
    % Assume the ID fo dA is the same as swc matrix 

    if size(root_node_ids) ~= size(parent_node_ids)
        error('The number of root node IDs must match the number of parent node IDs.');
    end

    %% STEP I: Identify floating components and node IDs in the compartment
    [node_compart_idx,num_components] = identify_connected_components(dA);
    if num_components ~= length (root_node_ids) + 1
        warning ('Number of floating branches different from input nodes pairs')
    end
    
    NodeID_old = swc_data(:,1); % original node ID 
    
if ~isempty (root_node_ids) % if is empty: reassign root node    
    if any (ismember ([root_node_ids, parent_node_ids ],NodeID_old)) 
        root_comp = node_compart_idx (root_node_ids);
        root_map = containers.Map(root_node_ids, root_comp);
        parent_comp = node_compart_idx (parent_node_ids);
        if all(sum(unique([root_comp', parent_comp'], 'rows') ~= [root_comp', parent_comp'], 2) == 2)
            error ('Child-parent connection within the same tree. Check input')
        end
    else
        error ('root or parent not within swc ID')
    end
else
    if  num_components == 1
        root_map = containers.Map(true_root_id, node_compart_idx(true_root_id));
    else
        error ('No re-assigining root node if the tree is not fully connected.')
    end
end

    if ~ismember (true_root_id,root_node_ids) 
        TrueRoot_comp = node_compart_idx (true_root_id);
        root_map (true_root_id) = TrueRoot_comp; % rootmap cotains all the root in all compartments
    else
        error ('The true root node should not be root node for a floating branch since it dosnt connect to anyone');
    end


    %% STEP 2: Sort floating branch so new root becomes assgined root (no longer lower-triag)
    dA_flipped = dA;
    allroot = [root_node_ids;true_root_id]; % also sort the compartemnet with the true root node
    for i = 1:length (allroot) 
        temp_root= allroot(i); % new node
        n_comp = root_map (temp_root);
        temp_node = NodeID_old(node_compart_idx == n_comp); % node ID in the floating branch        
        [dA_flipped] = changeRootAndFlipEdges(temp_node, dA_flipped, temp_root);
    end
    
    % Verify flipping
    temp = digraph (dA_flipped);
    outD = outdegree(temp);  % Get the out-degree of each node in the subgraph
    tempNode = NodeID_old( outD == 0); 
    if ~any(ismember (tempNode,allroot))
        error ('At least one of the floating branch is not flipped properly');
    else
        disp ('----- All floating branches flipped directions according to assigned root ---')
    end

    %% STEP 3: Directly connect child branch to other nodes (no need to comply to BCT formalism; sort later)
if ~isempty (root_node_ids) 
    pairs = [root_node_ids,parent_node_ids];
    dA_flipped_connected = dA_flipped;
    row = pairs(:,1); 
    col = pairs(:,2);
    for i = 1:length(row)
        dA_flipped_connected(row(i), col(i)) = 1;  % Assign each pair manually
    end 
    [~,num_comp_connected] = identify_connected_components(dA_flipped_connected);
    if num_comp_connected ~= 1
        error ('After connection all floating branhces, there are >1 comparement in the tree')
    end
else
    dA_flipped_connected = dA_flipped;
end

    %% STEP 4: Topologically sort the matrix
    G_connected = digraph(dA_flipped_connected'); % need to transpose the connected dA
    topoOrder = toposort(G_connected,"Order","stable");   
    dA_flipped_connected_sorted = dA_flipped_connected(topoOrder, topoOrder);
    figure; plot ( digraph(dA_flipped_connected_sorted),'Layout','layered');
    
    if isequal(triu(dA_flipped_connected_sorted, 1), zeros(size(dA_flipped_connected_sorted)))
        disp('--- The adjacency matrix is topolotically sorted and convert to lower triangular matrix --- .');
    else
        warning('The adjacency matrix is not lower triangular.');
    end

    % Convert dA to swc and visualize
    node_coords = swc_data(:,3:5);
    nodeR = swc_data(:,6);
    node_coords = node_coords(topoOrder,:);
    nodeR = nodeR (topoOrder);

    swc_data_updated = AdjMatrix2SWC (dA_flipped_connected_sorted,node_coords,nodeR);
    visualize_swc_interactive(swc_data_updated, 'NA');
end

%% Subfunctions
% Identify Connected Components
function [component_indices,num_components] = identify_connected_components(dA)
    G_undirected = digraph(dA + dA');  % undirected 
    component_indices = conncomp(G_undirected,'Type', 'strong');
    num_components = max(component_indices);
    fprintf('Number of connected components (including floating branches): %d\n', num_components);
    
    figure;
    G_directed = digraph (dA);
    plot ( G_directed,'Layout','layered');
    title ('Directd graph of the tree');
end

function swc_data = AdjMatrix2SWC (dA,node_coords,nodeR)
    % convert dA to swc table
    N = size(dA, 1);
    if size(dA, 2) ~= N
        error('Adjacency matrix dA must be square.');
    end
    
    swc_data = zeros(N, 7);
    swc_data(:, 1) = (1:N)';

    swc_data(:, 2) = 3;  % 3 = Default to dendrite type
    swc_data(:, 3:5) = node_coords;
    swc_data(:, 6) = nodeR;
    
    % Initialize parent IDs to -1
    swc_data(:, 7) = -1;  % Default parent ID is -1 (root)
    
    % Determine parent IDs from adjacency matrix
    [child_indices, parent_indices] = find(dA); % Root node no child --> sum (col) = 0 % Parentidx < child
    swc_data(child_indices, 7) = parent_indices;
end


function [dA_flipped] = changeRootAndFlipEdges(nodeIDs, dA, newRootID)
    % changeRootAndFlipEdges: flip the necessary edges for directed edges
    % after changing root for each floating branch
    
    figure; subplot (2,2,1);    
    G = digraph(dA ); 
    plot (G,'Layout','layered');

    subG = subgraph(G, nodeIDs);  % Create the subgraph of the selected nodes
    subplot (2,2,2);
    plot (subG,'Layout','layered');
    Map_Sub2Main = containers.Map(nodeIDs,1:numnodes(subG));

    % Find the current root node (node with out-degree 0 in the subgraph)
    outD = outdegree(subG); 
    oldRootIdx = find(outD == 0);  % Root has out-degree = 0
    oldRootID = nodeIDs(oldRootIdx);  % Convert index to actual node ID
    
    if isempty(oldRootIdx)
        error('No root node found with out-degree = 0 in the subgraph.');
    elseif length(oldRootIdx) > 1 
        warning('Multiple root nodes found. Taking the first one.')
    elseif min (nodeIDs) ~= oldRootID
        warning ('Old root ID is not the smallest of all the nodes ID in this floating branch for root %d',newRootID);
    end   

    subG_dA = adjacency (subG);
    UG = graph(subG_dA + subG_dA.');

    n = numnodes(UG);
    visited = false(n,1);
    pred = zeros(n,1);
    
    queue = Map_Sub2Main(newRootID);
    visited(Map_Sub2Main(newRootID)) = true;
    
    while ~isempty(queue)
        u = queue(1);
        queue(1) = [];
        neighbors_u = neighbors(UG, u);
        for i = 1:length(neighbors_u)
            v = neighbors_u(i);
            if ~visited(v)
                visited(v) = true;
                pred(v) = u;
                queue(end+1) = v; % pred = 0 means no parent
            end
        end
    end
    
    % Flip edges based on BFS traversal to point toward the new root
    G_sub_new = digraph(); % Initialize a new digraph
    G_sub_new = addnode(G_sub_new, n);
    
    for v = 1:n
        if pred(v) ~= 0 % If not the root
            u = pred(v); % Parent node
            % Edge from v (child) to u (parent) to point toward the root
            G_sub_new = addedge(G_sub_new, v, u);
        end
    end
    subplot (2,2,4);
    plot (G_sub_new,'Layout','layered'); 

    % map back the ID of subtree to main ID in dA
    sub_adjMatrix = adjacency(G_sub_new);
    sub_ID = cell2mat(values(Map_Sub2Main))';
    og_ID = cell2mat(keys(Map_Sub2Main))';
    [child, parent, ~] = find(sub_adjMatrix);
    newPair = zeros (length (child),2);
    for i = 1:length (child)
        c = og_ID(find(sub_ID== child(i), 1));
        p = og_ID(find(sub_ID== parent(i), 1));
        newPair (i,:) = [c,p];
    end

    dA_flipped = dA;
    dA_flipped (nodeIDs,:) = 0;
    for i = 1:size(newPair,1)
        dA_flipped(newPair(i,1), newPair(i,2)) = 1;  % Assign each pair manually
    end   
    
    G_flipped = digraph(dA_flipped ); 
    subplot (2,2,3);
    plot (G_flipped,'Layout','layered'); 
    pause (1); close gcf
end
