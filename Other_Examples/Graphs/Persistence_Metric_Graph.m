clear; close all; clc;

addpath('../../functions/persistence/');

% Save file info
save_to_file = false;
results_file = sprintf('Metric-Graph-%s', date());
image_file   = sprintf('Metric-Graph-%s', date());

% Parameters for the simulation
k = 1;      % Dimension of homology to compute
n = 2*k+2;  % Number of points in each configuration

% Number of n points samples to take
nReps = 1e5;

% List to save the resultsnReps
confs = zeros(n,3,nReps);
dms = zeros(n,n,nReps);
bd_times = zeros(nReps,2);

%% Create a graph by supplying the adjacency matrix
% Here are several examples

% % Two isoceles triangles whose bases are joined by two edges, each of
% % length 0.5. The bases of the triangles have length 1 and 4.
% A = [0,  0.5, 1,     0,   0, 0.5;
%     0.5, 0,   0.5,   0,   0,   0;
%     1,   0.5, 0,     0.5, 0,   0;
%     
%     0,   0,   0.5,   0,   1,   4;
%     0,   0,   0,     1,   0,   1;
%     0.5, 0,   0,     4,   1,   0];

% % Two triangles pasted at an edge
% A = [0, 1, 1, 0;
%      1, 0, 1, 1;
%      1, 1, 0, 1;
%      0, 1, 1, 0];

% % A cycle of length 2 (the cycle is formed by the vertices 1,2,3,4,5,6)
% % with an additional edge of length 0.9 between the vertices 2 and 5.
% A = [0,    0.4,  0,    0,   0,   0.5;
%      0.4,  0,    0.1,  0,   0.9, 0;
%      0,    0.1,  0,    0.5, 0,   0;
%      0,    0,    0.5,  0,   0.4, 0;
%      0,    0.9,  0,    0.4, 0,   0.1;
%      0.5,  0,    0,    0,   0.1, 0];

% % Two cycles of different lengths pasted at an edge
% A = [0, 1, 0, 0.5, 0, 0, 1;
%      1, 0, 1,   0, 0, 0, 0;
%      0, 1, 0,   1, 0, 0, 0;
%    0.5, 0, 1,   0, 1, 0, 0;
%      0, 0, 0,   1, 0, 1, 0;
%      0, 0, 0,   0, 1, 0, 1;
%      1, 0, 0,   0, 0, 1, 0];

% % Cycle of length 6 with 4 edges attached
% A = [0, 1, 0, 0, 0, 1, 1, 0, 0, 0;
%      1, 0, 1, 0, 0, 0, 0, 1, 0, 0;
%      0, 1, 0, 1, 0, 0, 0, 0, 0, 0;
%      0, 0, 1, 0, 1, 0, 0, 0, 1, 0;
%      0, 0, 0, 1, 0, 1, 0, 0, 0, 1;
%      1, 0, 0, 0, 1, 0, 0, 0, 0, 0;
%      1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
%      0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
%      0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
%      0, 0, 0, 0, 1, 0, 0, 0, 0, 0];

% % Cube
% A = [0, 1, 0, 1, 1, 0, 0, 0;
%      1, 0, 1, 0, 0, 1, 0, 0;
%      0, 1, 0, 1, 0, 0, 1, 0;
%      1, 0, 1, 0, 0, 0, 0, 1;
%      1, 0, 0, 0, 0, 1, 0, 1;
%      0, 1, 0, 0, 1, 0, 1, 0;
%      0, 0, 1, 0, 0, 1, 0, 1;
%      0, 0, 0, 1, 1, 0, 1, 0];

% Wedge of two cycles of different lengths
A = [0, 1.1, 1.1,   0, 0, 0;
     1.1, 0, 1.1,   0, 0, 0;
     1.1, 1.1, 0,   1, 0, 1;
     
     0, 0, 1,       0, 1, 0;
     0, 0, 0,       1, 0, 1;
     0, 0, 1,       0, 1, 0];

% % Wedge of two cycles of equal length
% A = [0,1,1,0,0;
%      1,0,1,0,0;
%      1,1,0,1,1;
%      0,0,1,0,1;
%      0,0,1,1,0];

% % Tree of cycles
% A = [0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0;
%      1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
%      0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
%      0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
%      0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
%      1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0;
%      0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0;
%      0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0;
%      0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0;
%      0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0;
%      1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1;
%      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1;
%      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0];

% MetricGraph defaults to a non-directed graph with uniform measure
% We can choose uniform or uniform_edges for the sampling method.
G = MetricGraph(A);

% See the graph
figure
G2 = graph(G.A);  % Convert my graph to a Matlab's graph object
p = plot(G2,'Layout','force');
% Uncomment the next line to have the edgelengths be proportional to their
% weight
% layout(p, 'force', 'WeightEffect','direct');

%% Execute the experiment
for i=1:nReps
    % -------------------------------------------------------------------------
    % Sample N points from the metric graph
    % -------------------------------------------------------------------------
    [X,dm] = sample(G, n);

    confs(:,:,i) = X;
    dms(:,:,i) = dm;

    % Find the persistence diagram
    [tb,td] = persistence_matrix(dm);
    bd_times(i,:)=[tb,td];
end % Of one of the nReps repetitions

% Save the results
if save_to_file
    save(sprintf('results/%s.mat',results_file)); %#ok<UNRCH>
end


%% See the results
if save_to_file
    f = figure('visible', 'off');
else
    f = figure('visible', 'on');
end

% Find the first Betti number to write a nice title
g = 1;
B1 = Betti_1(G);
name = sprintf('Graph %i, B1=%i, $\\mathbf{D}_{%i,%i}^{\\mathrm{VR}}(G)$', g, B1, n,k);
sgtitle(name, 'interpreter','latex');

% --------------------------------------
% Plot the graph using MATLAB's methods
% --------------------------------------
subplot(1,2,1);
G2 = graph(G.A);  % Convert my graph to a Matlab's graph object
p = plot(G2,'Layout','force');
layout(p, 'force', 'WeightEffect','direct');


% --------------------------------------
% Plot the persistence set
% --------------------------------------
subplot(1,2,2);
hold on;

% Main diagonal
tt=0:0.01:max(bd_times(:));
plot(tt,tt,'red');

% Plot results
plot(bd_times(:,1), bd_times(:,2), 'b.')
hold off;
