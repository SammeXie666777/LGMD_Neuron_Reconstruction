function export_hoc(swcdata, TypeName, outputname)
    % Perform export to NEURON format based on SWC file data.
    % Input:
    % swcdata - Nx7 matrix representing the SWC file
    % TypeName - cell array of section names corresponding to unique type IDs
    % outputname - string name for the output HOC file
    
    [~, cell_name, ~] = fileparts(outputname);
    fid = fopen(outputname, 'w');
    % ------------- Header
    % Error handling for file opening
    if fid == -1
        error('Could not open file for writing: %s', outputname);
    end

    % Write header information to the HOC file
    fprintf(fid, '/* --- Generated on: %s --- */\n\n', date());
    fprintf(fid, '/* ------------------------------------------\n');
    fprintf(fid, 'LGMD morphology swc file taken from Vaa3D Trace\n');
    fprintf(fid, 'Full reconstruction taken from 2-photon scans\n');
    fprintf(fid, 'Neuron name: %s\n', cell_name);
    fprintf(fid, '------------------------------------------*/\n\n');
    
    fprintf(fid, 'strdef neuron_name\n');
    fprintf(fid, 'neuron_name = "%s"\n\n', cell_name);
    
    fprintf(fid, '/* Create all sections */\n\n');
    
    % ------------- Create sections
    TypeID = sort(unique(swcdata(:, 2)),'ascend'); % Find unique types
    if numel(TypeName) ~= length(TypeID)
        error('Mismatch between number of section names and unique Type IDs');
    end

    % Ensure first node is soma
    if swcdata(1, 2) ~= 1 || swcdata(1, 7) ~= -1
        error('First node must be the soma and should have no parent');
    end   

    % Write sections to file
    for i = 1:numel(TypeName)
        %fprintf(fid, 'create %s[%d]\n',TypeName{i}, sum (swcdata(:,2) == TypeID (i)));
        numsegs = sum (swcdata(:,2) == TypeID (i));
        if TypeID (i) == 1
            numsegs = numsegs - 1;
        end
        fprintf(fid, 'create %s[%d]\n',TypeName{i}, numsegs);
    end

    % Create neuron index
    NEURON_Idx = cell(size(swcdata, 1), 2); 
    for i = 1:size(swcdata, 1)
        type = swcdata(i, 2);
        
        typeNameIdx = TypeID == type; 
        typeName = TypeName{typeNameIdx};
        
        % Calculate the NEURON index for this segment within its type
        typeSegmentIndices = find(swcdata(:, 2) == type);  % Get all indices of this type in swcdata
        segmentIndexWithinType = find(typeSegmentIndices == i) - 1;  % Find this segment's position within type
        NEURON_Idx{i, 1} = typeName;
        if i == 1 || i== 2
            NEURON_Idx{i, 2} = 0;
        else
            NEURON_Idx{i, 2} = segmentIndexWithinType;
        end
    end

    fprintf (fid, 'nseg = %d\n', 1);
    fprintf(fid, '\n');
    
    % ------------- Create geometry and connections
    fprintf(fid, '/* Geometry */\n');
    n_segs = length (swcdata(:,1));
    connections = cell(n_segs - 2 ,4); % n - 2 connections (first 2 are roots)
    
    % Create soma (SIZ)
    soma1X = swcdata(1, 3); soma2X = swcdata(2, 3);
    soma1Y = swcdata(1, 4); soma2Y = swcdata(2, 4);
    soma1Z = swcdata(1, 5); soma2Z = swcdata(2, 5);
    soma1R = swcdata(1, 6) * 2; soma2R = swcdata(1, 6) * 2;
    fprintf(fid, '\n');
    fprintf(fid, 'soma[%d] {\n', 0);
    
    fprintf(fid, '  pt3dclear()\n');
    fprintf(fid, '  pt3dadd(%5.2f, %5.2f, %5.2f, %5.2f)\n',soma1X, soma1Y, ...
    soma1Z, soma1R);
    fprintf(fid, '  pt3dadd(%5.2f, %5.2f, %5.2f, %5.2f)\n',soma2X, soma2Y, ...
    soma2Z, soma2R);
    fprintf(fid, '}\n');
    % Create all the other compartments
    for n = 3:size(swcdata, 1)   
        S2node_id = swcdata(n, 1); % current node ID 
        S2X = swcdata(S2node_id, 3);         % X coordinate
        S2Y = swcdata(S2node_id, 4);         % Y coordinate
        S2Z = swcdata(S2node_id, 5);         % Z coordinate
        S2R = swcdata(S2node_id, 6);         % Radius
        S2_NEURON_index = NEURON_Idx{S2node_id, 2};


        S1node_id = swcdata(S2node_id, 7);   % Parent node ID
        S1X = swcdata(S1node_id, 3);         % X coordinate
        S1Y = swcdata(S1node_id, 4);         % Y coordinate
        S1Z = swcdata(S1node_id, 5);         % Z coordinate
        S1R = swcdata(S1node_id, 6);         % Radius
        S1_NEURON_index = NEURON_Idx{S1node_id, 2}; % Neuron based index

        % Get the section name for this node type
        S2section_name = NEURON_Idx{S2node_id, 1}; % section name for child
        S1section_name = NEURON_Idx{S1node_id, 1}; % section name for parent
        % Write pt3dadd command for geometry: need two points - one child + one parent        
        fprintf(fid, [S2section_name '[%d] {\n'], S2_NEURON_index);
        fprintf(fid, '  pt3dclear()\n');
        fprintf(fid, '  pt3dadd(%5.2f, %5.2f, %5.2f, %5.2f)\n', S1X, S1Y, S1Z, 2 * S1R);  % Diameter = 2 * radius
        fprintf(fid, '  pt3dadd(%5.2f, %5.2f, %5.2f, %5.2f)\n', S2X, S2Y, S2Z, 2 * S2R);  % Diameter = 2 * radius       
        fprintf(fid, '}\n');

        % If the node has a parent, save the connection info
        connections{n-2, 1} = S2section_name;
        connections{n-2, 2} = S2_NEURON_index;
        connections{n-2, 3} = S1section_name;
        connections{n-2, 4} = S1_NEURON_index;
        
    end
    
    %--------------- Write the connections between sections
    fprintf(fid, '\n/* Connections */\n');
    for i = 1:size(connections, 1)
        child_label = connections{i, 1};
        child_index = connections{i, 2};
        parent_label = connections{i, 3};
        parent_index = connections{i, 4};
        
        if strcmp(child_label,'Axon') && strcmp(parent_label,'soma')
             choice = 0; % only for connection to soma
        else 
             choice = 1;
        end
        % Note: dend[0](0) refers to a point at the proximal end of dendrite 0, 
        % dend[0](1) refers to a point at the distal end of dendrite 0, and dend[0](0.2) refers to a point 20% down dendrite 0:       
        % hoc syntax: connect dend[0](0), soma[0](0)
        % fprintf(fid, 'connect %s[%d](%d), %s[%d](%d)\n',...
        %     parent_section, parent_index,choice, child_section, child_index,0);
        fprintf(fid, [parent_label '[%d] connect ' ...
                child_label '[%d](%f),%f\n'], ...
        parent_index, child_index, 0, choice);
    end
    
    % Close the file
    fclose(fid);
    
    fprintf('HOC file "%s" has been generated successfully.\n', outputname);
end
