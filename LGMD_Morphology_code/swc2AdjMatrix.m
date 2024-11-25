function [dA, node_coords, nodeR] = swc2AdjMatrix(swc_data)
    % Extract data from SWC table
    node_ids = swc_data(:, 1);
    x = swc_data(:, 3);
    y = swc_data(:, 4);
    z = swc_data(:, 5);
    nodeR = swc_data(:, 6);
    parent_ids = swc_data(:, 7);

    % Ensure node IDs are unique
    [unique_node_ids, ~, ~] = unique(node_ids);
    N = length(unique_node_ids);

    % Map node IDs to indices (1:N)
    [~, node_id_to_idx] = ismember(node_ids, unique_node_ids);
    [~, parent_id_to_idx] = ismember(parent_ids, unique_node_ids);

    % Identify valid parent-child relationships (excluding root nodes)
    valid_parents = parent_ids ~= -1;
    parent_indices = parent_id_to_idx(valid_parents);
    child_indices = node_id_to_idx(valid_parents);

    % Build the sparse adjacency matrix (directed from parent to child)
    dA = sparse(child_indices,parent_indices, 1, N, N);

    % Reorder node coordinates and radii to match unique_node_ids
    [sorted_node_ids, ~] = sort(unique_node_ids);
    [~, idx_in_node_ids] = ismember(sorted_node_ids, node_ids);

    node_coords = [x(idx_in_node_ids), y(idx_in_node_ids), z(idx_in_node_ids)];
    nodeR = nodeR(idx_in_node_ids);

    % Plot the adjacency matrix
    close all;
    spy(dA);
    title('Adjacency Matrix BCT format');
end
