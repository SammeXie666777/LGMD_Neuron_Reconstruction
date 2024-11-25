function [descendants] = FindNode_below(swcinfo,NodeStart,ExclusiveYN)
% Find node ID below some specified nodes topologically

    descendants = cell (1,length (NodeStart));

    for i = 1: length (NodeStart)
        top = NodeStart (i);
        descendants{i} = find_descendants(swcinfo, top);
        descendants{i} = [descendants{i};top];
    end

    % compare uniqueness of nodes from each parent 
    if strcmp(ExclusiveYN, "Y")
        % Ensure desendants are mutully exclusive
        is_unique = true; 
        for i = 1:numel(descendants)
            current_element = descendants{i};
            
            % Iterate over the rest of the cell elements
            for j = 1:numel(descendants)
                if i ~= j
                    other_element = descendants{j};
                    if isequal(current_element, other_element)
                        is_unique = false;
                        break;
                    end
                end
            end
            
            if ~is_unique
                break;
            end
        end
        
        if is_unique
            disp('All arrays in the cell are unique.');
        else
            disp('Some arrays in the cell are not unique.');
        end
    end
    
end

% Recursive function: Find all descendant nodes below the specified parent node.
function descendants = find_descendants(swc_data, parent_id)

    descendants = [];
    
    % Find the children of the specified parent_id
    children = swc_data(swc_data(:, 7) == parent_id, 1);

    % Add the children to the list of descendants
    descendants = [descendants; children];

    % Recursively find all descendants of each child
    for i = 1:length(children)
        child_descendants = find_descendants(swc_data, children(i));
        descendants = [descendants; child_descendants];
    end
end