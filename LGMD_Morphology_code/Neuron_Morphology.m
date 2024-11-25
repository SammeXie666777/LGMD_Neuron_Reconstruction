function [totalLength, meanRadius, surfaceArea] = Neuron_Morphology(swcData)
    % Input: Full path to the swc file containing the tree
    % Output: 3 parameters and results will be printed in console

    % read eswc
    %data = readmatrix(eswcFile, 'FileType', 'text', 'CommentStyle', '##');
    data = swcData;

    id = data(:, 1);
    x = data(:, 3);
    y = data(:, 4);
    z = data(:, 5);
    radius = data(:, 6);
    parent = data(:, 7);

    totalLength = 0;
    surfaceArea = 0;
    
    % Calculate the total length and surface area
    for i = 1:length(id)
        if parent(i) > 0
            % Find the parent node
            parentIndex = find(id == parent(i));
            if ~isempty(parentIndex)
                % Calculate distance between current node and its parent
                dx = x(i) - x(parentIndex);
                dy = y(i) - y(parentIndex);
                dz = z(i) - z(parentIndex);
                distance = sqrt(dx^2 + dy^2 + dz^2);
                
                % total distance
                totalLength = totalLength + distance;                
                % surface area: truncated cone
                r1 = radius(i);
                r2 = radius(parentIndex);
                surfaceArea = surfaceArea + pi * (r1 + r2) * sqrt(distance^2+(r1-r2)^2); 
            end
        end
    end
    % mean r
    meanRadius = mean(radius);

    fprintf('Total Length: %.2f(um)\n', totalLength);
    fprintf('Mean Radius: %.2f(um)\n', meanRadius);
    fprintf('Surface Area: %.2f(um^2)\n', surfaceArea);
end

