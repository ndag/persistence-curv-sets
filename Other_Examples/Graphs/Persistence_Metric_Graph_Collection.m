clear; close all; clc;

addpath('../../functions/persistence/');
addpath('../../functions/sampling/');

% In this experiment, we will sample nPoints points from each graph to create
% a test metric space X. From there, we will sample n points and calculate
% k-dim persistent homology for nReps times. We codify the sample with a
% index list I. We will save the test spaces in confs and dms, and include
% an extra list I_list to save the index lists for each graph.

% Name of the text file that contains the graphs
file_name = 'shortlist';
use_text_file = true;
adj_mats = read_graphs(file_name, use_text_file);

nGraphs = numel(adj_mats);
graphs = cell(nGraphs,1);

% Parameters for the simulation
k = 1;                  % Dimension of homology to compute
n = 2*k+2;              % Number of points in each configuration
nPoints = 100;          % Initial pool of points. From here we sample n points each time
nReps = 5000;           % Number of n point samples to search

% Save file info
save_to_file = false;

results_folder = sprintf('results/%s_%s', datestr(now, 'yy-mm-dd'), file_name);
if ~exist(results_folder, 'dir') && save_to_file
   mkdir(results_folder)
end

% Lists to save results --------------------------
confs = zeros(nPoints,3,nGraphs);
dms = zeros(nPoints,nPoints,nGraphs);
I_list = zeros(nReps,n,nGraphs);
bd_results = zeros(nReps,2,nGraphs);

% Repeat the experiment for each graph
for g=1:nGraphs
    % Instance a MetricGraph object using the adjacency matrix
    A = adj_mats{g};
    G = MetricGraph(A);
    graphs{g} = G;
    
    % Sample our test metric space X
    [X,dm_X] = sample(G, nPoints);
    
    confs(:,:,g) = X;
    dms(:,:,g) = dm_X;
    
    % Now we take random samples from X and find the persistent homology of
    % the sample
    bd_times = zeros(nReps,2);
    for i=1:nReps
        I = randi(nPoints,1,n);
        I_list(i,:,g) = I;
        
        % Find the distance matrix of this sample
        dm = dm_X(I,I);

        % Find the persistence diagram
        [tb,td] = persistence_matrix(dm);
        bd_times(i,:)=[tb,td];
    end % Of one of the nReps repetitions
    % Save the bd_times of this graph
    bd_results(:,:,g) = bd_times;
% end
% 
% %% Plotting
% for g=1:nGraphs
    A = adj_mats{g};
    G = MetricGraph(A);
    bd_times = bd_results(:,:,g);
    
    % See the results
    if save_to_file
        f = figure('visible', 'off');
    else
        f = figure('visible', 'on');
    end
    
    % Find the first Betti number to write a nice title
    B1 = Betti_1(G);
    name = sprintf('Graph %i, B1=%i', g, B1);
    sgtitle(name);
    
    % --------------------------------------
    % Plot the graph using MATLAB's methods
    % --------------------------------------
    subplot(1,2,1);
    G2 = graph(G.A);  % Convert my graph to a Matlab graph object
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
    
    % --------------------------------------
    % Save the results to a file
    % --------------------------------------
    if save_to_file
         x_width=8;  y_width=4;
        set(gcf, 'PaperPosition', [0 0 x_width y_width]);
        saveas(gcf, sprintf('%s/%s.png', results_folder, name));  %#ok<UNRCH>
    end
end % Of one of the nGraphs repetitions

% Save the results
if save_to_file
    save(sprintf("%s/data.mat",results_folder),...
        'adj_mats', 'bd_results', 'nPoints', 'confs', 'I_list', ...
        'file_name', 'n', 'k', 'nGraphs', 'nReps'...
    ); %#ok<UNRCH>
end
