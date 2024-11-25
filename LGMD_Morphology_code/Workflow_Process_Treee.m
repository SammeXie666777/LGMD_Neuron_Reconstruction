%% Implementing visualization & morphology calcualtion & Edit connection & change segment types & Export swc
% Make sure to add the following functions to your path:
%% 
% * visualize_swc_interactive.m
% * Neuron_Morphology.m
% * process_floating_branches.m

% Get swc file
[swcfilename, pathname] = uigetfile({'*.swc;*.eswc;*.mat', 'SWC and ESWC Files (*.swc, *.eswc)'});
swcpath = fullfile(pathname,swcfilename);
[~,~,ext] = fileparts (swcpath);
if contains (ext,'mat') % if the data is in the format of .mat
    load (swcpath,"swc_data_updated"); % make sure only one variable containing swc marix 
else                    % if the data is in the format of .swc (text file)
    swc_data_updated = readmatrix(swcpath, 'FileType', 'text', 'CommentStyle', '##');
end

%%%%%%%%%%%%%%%%%%%
%% Visualization and morphology
% Visualization
visualize_swc_interactive(swc_data_updated, 'NA'); 
% Morphology calculation
[Length, MeanR, surfaceArea] = Neuron_Morphology (swc_data_updated);

%%%%%%%%%%%%%%%%%%
%% Edit connection
% Delete previous connection and connect to new nodes to make a fully connected 
% tree
% Step 1:  Use TREE toolbox to load the tree and visualize current Tree topology

    trees = load_tree(swcpath); % calls repair_tree and sort_tree according to BCT formalism
    dA_tree (trees);
    title ('Adjacency matrix dA for sorted, original tree. x = parent; y = child; X<Y');
    dA = trees.dA;
    node_coords = [trees.X, trees.Y, trees.Z];
    nodeR = trees.D/2;
    swc_data = AdjMatrix2SWC (dA,node_coords,nodeR);
    disp('Original tree; Highlight nodes for future reference')

    visualize_swc_interactive(swc_data, 'Tip')

    % G_undirected = digraph(dA );  % Symmetrize the adjacency matrix
    % plot ( G_undirected,'Layout','layered');
    % component_indices = conncomp(G_undirected);
    % num_components = max(component_indices);
    % fprintf('Number of connected components (including floating branches): %d\n', num_components);

%% Step 2: Input node id_A that disconnected from id_B and connect to id_A to id_C


    % Potentially write a GUI that allows me to 1) select pairs using
    % different color 2) select pairs directly delete the connections 
    % PairtoOld = [1494,1495; 1467,1466];

    pairIDold = [1503,1478; 1474,1473; 1471,1470; 41,42; 27,28; 78,79; 212,211; 120,110; 705,706; 
        230,231;230,229; 176, 177; 128, 129; 1235, 89; 824,823; 120, 970; 1000,1001;
        984,983; 2617, 2616; 985,984; 969,120;51,52; 539,199; 1161,1162; 450,451; 574,573;
        558,557; 574,573]; % Pairs of ID to be disconnedted 
    pairIDold = sort(pairIDold,2,"descend");
    xcoor = pairIDold(:,1); 
    ycoor = pairIDold (:,2);
    dA_mod = dA;
    dA_mod (xcoor,ycoor) = 0; % Diconnect
    swc_datav1 = AdjMatrix2SWC (dA_mod,node_coords,nodeR); % The function is at the bottom
    figure;
    spy(dA_mod); 
    visualize_swc_interactive(swc_datav1, 'NA') % visualize if it connections cut 

    
%% Step 3: Connect floating branches (not sorted and will result in trifurcation)

    % dA is below diagnol = parent id > child
    % root_node: new root for every floating branch
    % connect_to_node: node that each floating branch coonect to 
    % Sum (dA (i,:)) = num parents ([0,1]); Sum (dA (:,j)) = num children ([0,2])
    

    [filename, pathname] = uigetfile({'.fig'});
    disp ('Highlight nodes that comes the new root for each floating branch');
    [root_node_ids] = visualize_swc_interactive(swc_datav1, 'Tip',  fullfile(pathname,filename)); % Find the root node that will be the root for each floating branch

    %%
    [connectedTo] = visualize_swc_interactive(swc_datav1, 'Tip',  fullfile(pathname,filename),'r');

    %% Optional: make sure the true root node is not part of the 
    floatingRoot = root_node_ids(:,1);
    floatingRoot (floatingRoot == 3) = [];
    ConnectedNode = connectedTo (:,1); % exclude the true rooooot node for the entire tree from the ConnectedNode variables

%% Step 4
    [dA_flipped_connected_sorted, swc_data_updated] = process_floating_branches(dA_mod, swc_datav1, 1, floatingRoot, ConnectedNode);
    %save ('SWC_dA_connected.mat',"swc_data_updated","dA_flipped_connected_sorted");

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assign segment type to each node
% Sort the tree so the SIZ has ID 1 and 2: make the root of the tree 1

true_root_id = 2897;
dA = swc2AdjMatrix (swc_data_updated);
[~, swc_data_sorted] = process_floating_branches(dA, swc_data_updated, true_root_id, [], []);

swc_data_sorted (:,2) = 0;

%% Need to SORT FIRST
% Assign |*SIZ or soma = 1*|: 2 nodes define a compartment
disp ('Assign SIZ (2 nodes define a compartment)...')
SIZ_ID = [1,2];
swc_data_sorted (SIZ_ID,2) = 1;

%% Need the number of nodes within each type to be paired
% Assign |*Field A = 3*|, |*B = 4*|, |*C = 5*|, |* CellBody (real soma) = 6*|, |*D = 7*|
disp ('Assign Field A, B, C, D, and soma...')
NodeA_ID = 47;
NodeB_ID = 1906;
NodeC_ID = 2579;
NodeD_ID = 1819; % Arbitutary
CellBodyID = 1750;
AllTopNode = [NodeA_ID, NodeB_ID, NodeC_ID, CellBodyID, NodeD_ID]; % [2,3,4,5,6] 
[descendants] = FindNode_below(swc_data_sorted, AllTopNode, "Y");

% Assign Seg Type to each node
AssignedID = [3,4,5,6,7];
for i = 1:numel (descendants)
    swc_data_sorted(descendants{i},2) = AssignedID (i);
end

%% 
% Assign |*Axon = 2*|
idx = swc_data_sorted(:,2) == 0;
swc_data_sorted (idx,2) = 2;

disp ('Assigned all nodes... Adjust node types for indiviual nodes accordingly');
swc_data_sorted (1749,2) = 4;

ID = 1:7;
TypeName = {'soma','Axon','FieldA','FieldB','FieldC','Cellbody','FieldD'}; %[1 - 7]
for i = 1:numel (unique (swc_data_sorted(:,2)))
    fprintf ('Number of nodes for node type %s with ID %d n = %d\n',TypeName{i},  ID (i), sum (swc_data_sorted(:,2) == ID (i)));
end
%% Turn swc matrix to swc file
swcfilename = 'Fully_Connected_LGMD_sorted.swc';
save_to_swc(swc_data_sorted, swcfilename)

%% Turn the swc file into hoc file directly
TypeName = {'soma','Axon','FieldA','FieldB','FieldC','Cellbody','FieldD'}; %[1 - 7]
output_hoc = 'LGMD_Complete_Construction.hoc';
export_hoc(swc_data_sorted, TypeName, output_hoc)


%% Subfunction
% convert swc to ad matrix
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