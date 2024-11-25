function save_to_swc(data, filename)
    % Save neuron data to an SWC file
    % Input:
    %   data: Nx7 matrix where each row is [NodeID, Type, X, Y, Z, Radius, ParentID]
    %   filename: Name of the SWC file to write (e.g., 'neuron.swc')
    
    % Open file for writing
    fileID = fopen(filename, 'w');
    
    if fileID == -1
        error('Could not open file for writing.');
    end
    
    % Write SWC header (optional but often included for documentation)
    fprintf(fileID, '# SWC format generated from MATLAB on %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fileID, '# Columns: ID Type X Y Z Radius ParentID\n');
    
    % Write the data (one line per row of the matrix)
    for i = 1:size(data, 1)
        fprintf(fileID, '%d %d %.3f %.3f %.3f %.1f %d\n', ...
            data(i, 1), ...   % Node ID
            data(i, 2), ...   % Node type
            data(i, 3), ...   % X-coordinate
            data(i, 4), ...   % Y-coordinate
            data(i, 5), ...   % Z-coordinate
            data(i, 6), ...   % Radius
            data(i, 7));      % Parent ID
    end
    
    % Close the file
    fclose(fileID);
    
    disp(['SWC file saved as ', filename]);
end
