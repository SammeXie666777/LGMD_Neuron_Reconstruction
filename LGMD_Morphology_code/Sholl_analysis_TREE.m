start_trees; 
%% Load swc file
[filename, pathname] = uigetfile({'*.swc;*.eswc', 'SWC and ESWC Files (*.swc, *.eswc)'});
swcpath = fullfile(pathname,filename);
LGMD = load_tree(swcpath); % calls repair_tree
trees{1}.frustum = 1;
LGMD.frustum = 1;
%% Visual exploration
% xplore_tree; axis off;
figure;
plot_tree(LGMD,[1 0 0]);
shine; axis off; view(15,55);

% figure;
% xdend_tree (LGMD,'-s')

%% Stats/Morph
LGMDStats = stats_tree   (LGMD, [], [], '-x,-s'); % Global stats
seg_surf = surf_tree (LGMD);    % surface area
totalSurfA = sum(seg_surf);
seg_vol = vol_tree (LGMD);      % Volumn (sth wrong w )
totalVol = sum(seg_vol);
numTips = sum (T_tree (LGMD));
numBranch = sum (B_tree (LGMD))*2+1;

[spanning, ~] = gscale_tree(trees); % Spanning field of trees

% write stats in text file
fileID = fopen(fullfile(pathname,'LGMDStats_summary.txt'), 'w');
fprintf(fileID, 'Total cable length: %2d(um)\n', LGMDStats.gstats.len);
fprintf(fileID, 'Total surface area: %2d (um^2)\n', totalSurfA);
fprintf(fileID, 'Total volume: %2d (um^3)\n', totalVol);
fprintf(fileID, 'Average diameter: %2d (um)\n', mean(LGMD.D));
fprintf(fileID, 'Num of branches: %d\n', numBranch);
fprintf(fileID, 'Num of branch points: %d\n', LGMDStats.gstats.bpoints);
fprintf(fileID, 'Num of termination points: %d\n', numTips);

fclose(fileID);

%% Sholl
disp('Sholl analysis using sholl_tree function');
stepSize = 50; % increasing diameter for concentric sphere
[s, dd, sd, XP, YP, ZP, iD] = sholl_tree (LGMD,stepSize,'-3s'); % 3D
sholl_tree (LGMD,stepSize,'-s'); % 2D graph
% s: num intersection at diameters dd
% sd: double intersection
% XP, YP, ZP: coordinates of intersection pts
% iD: index of these points in dd
 
figure; % Plot num intersection wrt dd
plot(dd, s,'LineWidth',3,'Marker','o');
xlabel('Diameter (µm)');
ylabel('Num of intersections');
title(['Sholl Analysis with Step Size ', num2str(stepSize), ' µm']);

%% GUI
cgui_tree

%% tree topology for a graphical perspective
dA_tree (LGMD)

