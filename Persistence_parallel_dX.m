clear; close all; clc;
addpath('functions/persistence');

% Record the starting time
tic
% -------------------------------------------------------------------------
% Bookkepping variables
% -------------------------------------------------------------------------
CoresRequested = 0;         % Number of parallel workers used
% Set CoresRequested = 0 to use Matlab's default choice

save_to_file = true;
results_file = 'Persistence_parsave';

temp_file = 'parsave';

% -------------------------------------------------------------------------
% Parameters for the simulation
% -------------------------------------------------------------------------
nReps = 1e4;        % Number of snapshots
dim = 1;            % Dimension of persistence VR diagramns 
n = 2*dim+2;        % Number of points in each configuration

% Here load your distance matrix dX representing the metric
% space. Keep in mind that a reasonable size for dX is at most 5k x
% 5k (depending on your system).
filename_metric_space = 'dX_example_circle.mat';
load(filename_metric_space,'dX');

% -------------------------------------------------------------------------
% Set up the parallel pools and variables
% -------------------------------------------------------------------------
% Initialize the parallel pool if it doesn't exist yet
if isempty(gcp('nocreate'))
    if CoresRequested ~= 0
        parpool('local',CoresRequested);
    else
        parpool('local');
    end
end
p = gcp('nocreate');
nCores = p.NumWorkers;

% Create the directory temp if it doesn't exist
if ~exist('temp','dir')
    mkdir('temp');
end

% Global result variables
confs = zeros(n,nReps);
dms = zeros(n,n,nReps);
bd_times = zeros(nReps,2);

% Create a parallel.pool.Constant variable. The workers will use them to
% save results while the program runs. In case the program crashes, these
% files can be used to read intermediate results
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

% Change dX to a Constant object, so that all workers can access it
dX_constant = parallel.pool.Constant(dX);

t_setup = toc;   % Record the ending time of the setup
tic         % Record when the simulation itself started

% -------------------------------------------------------------------------
% Execute the experiment
% -------------------------------------------------------------------------
parfor i=1:nReps
    % Sample n points from dX
    dX_0 = dX_constant.Value;
    I = randi(length(dX_0),1,n);      % get n indices
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

    % Save results
    confs(:,i) = I';
    dms(:,:,i) = dm;
    bd_times(i,:) = [tb,td];
    
    % And save them to disk (for redundancy)
    matfileObj = myMatfileConstant.Value;
    matfileObj.confs(:,i) = I';
    matfileObj.dms(:,:,i) = dm;
    matfileObj.bd_times(i,:) = [tb,td];
    matfileObj.gotResult(1,i) = true;

    % Show progress
    if rem(i,floor(nReps*0.05))==0
        % progress = sum(matfileObj.gotResult,2);
        % fprintf(1,'Progress (%02d) %d/%d\n',t.ID, progress, nReps);
        t = getCurrentTask();
        progress = 100*sum(matfileObj.gotResult,2)/nReps;
        fprintf(1,'Progress %02.0d/%02.0d\n',i,nReps);
        fprintf(1,'Progress (%02d) %02.0d%%\n',t.ID, progress);
    end
end % Of one of the nReps repetitions

% Record the ending time of the simulation and display the duration of the
% experiment
t_calculations = toc;
fprintf("Preparations: %03f\n",t_setup);
fprintf("Calculations: %0.3f\n",t_calculations);

% -------------------------------------------------------------------------
% Graph results
% -------------------------------------------------------------------------
figure();
hold on;

tt = 0:0.01:max(bd_times);
plot(tt, tt, 'red');                        % Draw the diagonal for reference

plot(bd_times(:,1), bd_times(:,2), 'b.');   % Plot the persistence set
hold off;

% -------------------------------------------------------------------------
% Save results
% -------------------------------------------------------------------------
if save_to_file
    exclude = "t|p|myFname|myMatfile|myMatfileConstant|dX_constant|exclude";
    save('Persistence_ripser_parallel', '-regexp', sprintf('^(?!(%s)$).', exclude)); %#ok<UNRCH>
end

% Delete the parallel pool
% delete(gcp('nocreate'))