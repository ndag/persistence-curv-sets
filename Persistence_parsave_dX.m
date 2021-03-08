clear; close all; clc;

% Start a timer
tic

addpath(char("./functions"));

% set the flag below:
save_to_file = false;

results_file = 'Persistence_parsave';
temp_file = 'parallel';

% Parameters for the simulation
dim = 1;        % Dimension of persistence VR diagramns 
n = 2*dim+2;    % Number of points in each configuration

nCores = 6;         % Number of parallel workers used

nReps = 100;        % Number of snapshots


% Here load your distance matrix dX representing the metric
% space. Keep in mind that a reasonable size for dX is at most 5k x
% 5k (depending on your system).
filename_metric_space = './dX.mat';
load(filename_metric_space,'dX');

% Initialize the parallel pool
parpool('local',nCores);

t1 = toc;   % End timer for setup
tic         % Start a timer for the simulation

% Create a parallel.pool.Constant from the 'Composite'
% This allows the worker-local variable to be used inside PARFOR
spmd
    t = getCurrentTask();
    myFname = sprintf("%s/temp/checkpoint_%s_%02d", pwd, temp_file, t.ID); % each worker gets a unique filename
    myMatfile = matfile(myFname, 'Writable', true);
    
    % Seed the variables in the matfile object
    myMatfile.confs = zeros(n,nReps);
    myMatfile.dms = zeros(n,n,nReps);
    myMatfile.bd_times = zeros(nReps,2);
    myMatfile.gotResult = false(1,nReps);
end
myMatfileConstant = parallel.pool.Constant(myMatfile);

% Define dX as a Constant object, so that all workers can access it
dX = parallel.pool.Constant(dX);

% -------------------------------------------------------------------------
% Repeat the experiment
% -------------------------------------------------------------------------
t1 = toc;
tic
parfor i=1:nReps
    % -------------------------------------------------------------------------
    % Sample n points from dX
    % -------------------------------------------------------------------------
    dX_0 = dX.Value;
    
    I = randi(length(dX_0),n,1);         % get n indices
                                       % sampled from 1:length(dX)
    
    % Pick dm as a submatrix of dX that corresponds to the indices chosen
    % in I. This creates a metric space with the distance matrix dX
    % restricted to the n points in I
    P = nchoosek(I,2);
    inds = sub2ind([length(dX_0),length(dX_0)],P(:,1),P(:,2));
    dm = dX_0(inds);
    dm = squareform(dm);

    % Compute the Vietoris-Rips persistence diagrams of the metric space
    % with distance matrix dm
    [tb,td] = persistence_matrix(dm);
    bd_times(i,:) = [tb,td];

    % Save the results of the iteration
    t = getCurrentTask();
    matfileObj = myMatfileConstant.Value;
    
    matfileObj.confs(:,i) = I;
    matfileObj.dms(:,:,i) = dm;
    matfileObj.bd_times(i,:) = bd_times(i,:);
    matfileObj.gotResult(1,i) = true;

    % Show progress
    if rem(i,floor(nReps*0.05))==0
        % progress = sum(matfileObj.gotResult,2);
        % fprintf(1,'Progress (%02d) %d/%d\n',t.ID, progress, nReps);
        progress = 100*sum(matfileObj.gotResult,2)/nReps;
        fprintf(1,'Progress %02.0d/%02.0d\n',i,nReps);
        fprintf(1,'Progress (%02d) %02.0d%%\n',t.ID, progress);
    end
end % Of one of the nReps repetitions

t = toc;
fprintf("Preparations: %03f\n",t1);
fprintf("Calculations: %0.3f\n",t);

% Save the results
if save_to_file
    exclude = "myFname|myMatfile|myMatfileConstant|exclude";
    save(results_file, '-regexp', sprintf('^(?!(%s)$).', exclude));
end

% Delete the parallel pool
delete(gcp('nocreate'))