%% Visualization and selecting interested node
function [selected_points] = visualize_swc_interactive(data, TiporNA, figurePath, color)
    % Visualizes SWC data and allows interactive selection of nodes.
    %
    %   [selected_points] = visualize_swc_interactive(data, BrushorTiporNA, figurePath, color)
    %
    % Parameters:
    %   data           - SWC data matrix [index, type, x, y, z, radius, parent]
    %   TiporNA - Mode of selection: 'tip', or 'NA'
    %   figurePath     - Path to an existing figure file to open (optional)
    %   color          - Color to use for highlighting selected nodes (optional, default is 'g')
    %
    % Returns:
    %   selected_points - Matrix of selected nodes (SWC format)
    %                   - Saved as a .mat file with timestampled filename
    % Example call:
    %   [selected_points] = visualize_swc_interactive(SWC data matrix, 'tip', figurePath, color)

    % Set default color if not provided
    if nargin < 4 || isempty(color)
        color = 'g';  % Default color is green
    end

    % Check if figure path is provided and valid
    if nargin >= 3 && ~isempty(figurePath) && (ischar(figurePath) || isstring(figurePath))
        % Open the figure from the given file path
        fig = openfig(figurePath, 'reuse', 'visible');
        figure(fig);  % Set the opened figure as the current active figure
        disp(['Figure "', figurePath, '" is now the active figure.']);

        % Attempt to retrieve scatter plot data from the figure
        ax = gca;  % Get current axes
        scatter_plot = findobj(ax, 'Type', 'Scatter');
        if isempty(scatter_plot)
            error('No scatter plot found in the loaded figure.');
        end

        if nargin < 1 || isempty(data)  % If data input empty, try retrieving it from the figure's UserData
            data = get(fig, 'UserData');  % Retrieve stored SWC data
            if isempty(data)
                error('No SWC data found in the figure or as input.');
            end
        elseif ~ismember(length(data(1,:)), [7, 12])
            error('Wrong SWC data input format');
        end
    else
        % No figure path provided, plot the SWC data as a new figure
        fig = figure;
        ax = axes(fig);
        hold(ax, 'on');

        % Extract SWC data fields
        index = data(:, 1);
        x = data(:, 3);
        y = data(:, 4);
        z = data(:, 5);
        radius = data(:, 6);
        parent = data(:, 7);

        % Plot each segment of the neuron
        for i = 1:length(index)
            if parent(i) ~= -1
                parentIdx = find(index == parent(i)); % Parent index
                if ~isempty(parentIdx)
                    plot3(ax, [x(i), x(parentIdx)], [y(i), y(parentIdx)], [z(i), z(parentIdx)], 'k-', 'LineWidth', 2);
                else
                    warning('Parent index not found for node ID: %d', index(i));
                end
            end
        end

        % Plot the nodes
        scatter_plot = scatter3(ax, x, y, z, radius * 10, 'filled', 'MarkerFaceColor', 'auto', 'MarkerEdgeColor', 'k');

        % Plot the root node
        rootIdx = find(parent == -1);
        if ~isempty(rootIdx)
            scatter3(ax, x(rootIdx), y(rootIdx), z(rootIdx), radius(rootIdx) * 10, 'filled', 'MarkerFaceColor', 'blue');
        else
            warning('No root found in the tree');
        end

        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        title('Vaa3D Tracing Visualization');
        grid on;
        rotate3d on;  % Allow interactive rotation
        hold off;

        % Optionally store the data in the figure for future use
        set(fig, 'UserData', data);
    end

    % Initialize selected_points
    selected_points = [];

    % Determine the data type of node IDs
    keyType = class(data(1, 1));

    % Initialize highlightedNodes with the correct KeyType
    highlightedNodes = containers.Map('KeyType', keyType, 'ValueType', 'any');

    % Check whether the user chose 'tip' or 'brush' mode
    switch TiporNA
        case 'Tip'
            % Set the custom CloseRequestFcn
            set(fig, 'CloseRequestFcn', {@save_and_close_figure, fig});

            % Data Cursor Mode (Tip)
            disp('Select nodes through the cursor. Hit Enter to save current node and close figure to finish selection.');
            dcmObject = datacursormode(fig);
            set(dcmObject, 'DisplayStyle', 'datatip', 'SnapToDataVertex', 'on');
            set(dcmObject, 'UpdateFcn', {@customDataTipFcn, data(:, 3:5)});  % Custom cursor tip to display node ID
            datacursormode on;

            keep_running = true;
            while keep_running && ishandle(fig)
                pause;
                if ishandle(fig)
                    cursor = getCursorInfo(dcmObject);
                    if ~isempty(cursor)
                        newStat = [cursor.Position(1), cursor.Position(2), cursor.Position(3)];
                        [~, node_id] = min(vecnorm(data(:, 3:5) - newStat, 2, 2));  % Find node index

                        hold(ax, 'on');  % Ensure hold is on

                        % Get the key for the current node
                        key = data(node_id, 1);

                        if isempty(selected_points) || isempty(find(selected_points(:, 1) == key, 1))
                            % Highlight the selected node using the specified color
                            h = scatter3(ax, data(node_id, 3), data(node_id, 4), data(node_id, 5), 100, color, 'filled');
                            highlightedNodes(key) = h;
                            selected_points = [selected_points; data(node_id, :)];
                            fprintf('Node ID: %d, Coordinates: (%.2f, %.2f, %.2f)\n', key, data(node_id, 3), data(node_id, 4), data(node_id, 5));
                        else
                            % Un-highlight the node
                            if isKey(highlightedNodes, key)
                                delete(highlightedNodes(key));
                                remove(highlightedNodes, key);
                            else
                                warning('Key %d not found in highlightedNodes.', key);
                            end
                            fprintf('Node ID: %d is un-selected.\n', key);
                            idx = find(selected_points(:, 1) == key, 1);
                            selected_points(idx, :) = [];
                        end

                        % Clear the cursor info to avoid repeated selection
                        datacursormode off;
                        datacursormode on;
                    end
                else
                    keep_running = false;
                end
            end
            % Save selected_points at the end with timestamp
            if ~isempty (selected_points)
                timestamp = datestr(now, 'yyyymmdd_HHMMSS');
                filename = sprintf('saved_node_info_%s.mat', timestamp);
                save(filename, 'selected_points');
                disp(['Updated selected_points saved to ', filename]);
            end

        case 'NA'
            disp('Displaying tree image and nodeID cursor mode on');
            set(fig, 'CloseRequestFcn', 'closereq');
            
            % Enable Data Cursor Mode to display node IDs
            dcmObject = datacursormode(fig);
            set(dcmObject, 'DisplayStyle', 'datatip', 'SnapToDataVertex', 'on');
            set(dcmObject, 'UpdateFcn', {@customDataTipFcn, data(:, 3:5)});  % Custom cursor tip to display node ID
            datacursormode on;

        otherwise
            warning('Invalid selection mode specified. No interactive selection will occur.');
            set(fig, 'CloseRequestFcn', 'closereq');
    end

    %% --------- Nested Functions ---------
    % Custom Data Tip Function to Display Node ID
    function txt = customDataTipFcn(~, event_obj, myarray)
        pos = get(event_obj, 'Position');
        [~, node_id] = min(vecnorm(myarray - pos, 2, 2));
        txt = {['Node ID: ', num2str(data(node_id, 1))]};
    end

    function save_and_close_figure(~, ~, fig)
        % Prompt the user whether to save the figure
        choice = questdlg('Do you want to save the figure?', 'Save Figure', 'Yes', 'No', 'Yes');
        
        if isempty(choice)
            % User closed the dialog without making a selection
            disp('No selection made. Returning to the figure.');
            return;  % Do not proceed further; return control to the figure
        end
    
        switch choice
            case 'Yes'
                % Remove the custom CloseRequestFcn to prevent it from persisting
                set(fig, 'CloseRequestFcn', 'closereq');
                
                % Prompt the user for the save location
                [file, path] = uiputfile('*.fig', 'Save Figure As');
                if isequal(file, 0) || isequal(path, 0)
                    disp('User canceled the save operation.');
                else
                    % Save the figure to the specified location
                    try
                        savefig(fig, fullfile(path, file));
                        disp(['Figure saved as ', fullfile(path, file)]);
                    catch ME
                        warning(ME.identifier, '%s', ME.message);
                    end
                end
                % Close the figure
                delete(fig);
            case 'No'
                disp('Figure not saved.');
                % Remove the custom CloseRequestFcn to prevent it from persisting
                set(fig, 'CloseRequestFcn', 'closereq');
                % Close the figure
                delete(fig);
        end
    end

end
