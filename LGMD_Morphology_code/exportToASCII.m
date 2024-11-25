function exportToASCII(data, filename, label)
    % exportToASCII: Exports a two-column matrix to an ASCII file with a specific format
    %
    % Inputs:
    %   data     - A matrix with two columns: first column is x, second column is y
    %   filename - Name of the output ASCII file (e.g., 'output.txt')
    %   label    - Optional label string (e.g., 'soma.v(.5)')
    %
    % Example:
    %   data = [0, -70; 0.025, -65.8269; ...];
    %   exportToASCII(data, 'output.txt', 'LGMD')

    % Open the file for writing
    fid = fopen(filename, 'w');
    if fid == -1
        error('Cannot open file: %s', filename);
    end

    % Write the label line if provided
    if nargin > 2 && ~isempty(label)
        fprintf(fid, 'label:%s\n', label);
    end

    % Write the count line
    count = size(data, 1);
    fprintf(fid, '%d\n', count);

    % Write the data: one tab-separated x y data pair per line
    for i = 1:count
        fprintf(fid, '%f\t%f\n', data(i, 1), data(i, 2));
    end

    % Close the file
    fclose(fid);
end
