clear; close all; clc;

tic
save_to_file = false;

addpath(char("./funciones"));

results_file = 'Persistence_parsave';
temp_file = 'parallel';

% Parameters for the simulation
n = 6;      % Number of points in each configuration
d = 2;      % Dimension of the sphere
nCores = 2; % Number of parallel workers used

% Bump this up
nReps = 200;

% Parameters for the Rips complex
max_dimension = d;

% List to save the results
confs = zeros(n,d+1,nReps);
dms = zeros(n,n,nReps);
bd_times = zeros(nReps,2);

%Create a mat-file per worker using SPMD
% Initialize the parallel pool
% parpool('local',nCores)
t1 = toc;
tic
spmd
    t = getCurrentTask();
    myFname = sprintf("%s/temp/checkpoint_%s_%02d", pwd, temp_file, t.ID); % each worker gets a unique filename
    myMatfile = matfile(myFname, 'Writable', true);
    % Seed the variables in the matfile object
    myMatfile.confs = zeros(n,d+1,nReps);
    myMatfile.dms = zeros(n,n,nReps);
    myMatfile.bd_times = zeros(nReps,2);
    myMatfile.gotResult = false(1,nReps);
end
% Create a parallel.pool.Constant from the 'Composite'
% This allows the worker-local variable to be used inside PARFOR
myMatfileConstant = parallel.pool.Constant(myMatfile);

% -------------------------------------------------------------------------
% Repeat the experiment
% -------------------------------------------------------------------------
t1 = toc;
tic
parfor i=1:nReps
    % -------------------------------------------------------------------------
    % Sample n points from S^d
    % -------------------------------------------------------------------------
    X = randn(n,d+1);               % get n points in R^{d+1}
    nor = sum(X.^2,2).^0.5;         % get norm
    X = X./(nor*ones(1,d+1));       % normalize
    
    % Get the distance matrix
    dm = L2_distance(X',X');
    dm = 2*asin(dm/2);

    confs(:,:,i) = X;
    dms(:,:,i) = dm;

    % Create a metric space with the distance matrix of the n points in X
    % and compute its Vietoris-Rips complex
    t = getCurrentTask();
    PDs = RipsFiltrationDMPar(dm, max_dimension, t.ID, temp_file);

    % Get birth and death times
    % We only expect one interval in dimension 2, so we just get the first
    H_d = PDs{d+1};
    if numel(H_d) ~= 0
        % Add them to the list
        bd_times(i,:) = H_d(1,:);
    else
        % If there is no persistence, we record it
        bd_times(i,:) = [0,0];
    end
    
    % Save checkpoint
    matfileObj = myMatfileConstant.Value;
    
    matfileObj.confs(:,:,i) = X;
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
% delete(gcp('nocreate'))
