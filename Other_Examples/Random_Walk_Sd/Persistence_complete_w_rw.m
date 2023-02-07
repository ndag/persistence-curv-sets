% This script calculates a first approximation to a diagram D_n,k(S^d).
% Because we have seen that uniform sampling doesn't cover the whole
% diagram, we do a Random Walk after that to try to fix that.
clear; close all; clc;
tic

% Load functions and previous configurations
addpath('../../functions/persistence/');
addpath('../../functions/functions_S2/');

% -------------------------------------------------------------------------
% Parameters for the simulation
% -------------------------------------------------------------------------
k = 2;              % Dimension of homology to calculate
n = 2*k+2;          % Number of points in each configuration
d = 2;              % Dimension of the sphere
sigma = 0.075;      % Standard deviation of the normal perturbation
eps = 0.005;        % Radius of density neighborhood
rng(1);             % set random seed

% Number of configurations to sample
nReps = 5e3;
nSteps = 5e3;

% -------------------------------------------------------------------------
% Bookkepping variables
% -------------------------------------------------------------------------
CoresRequested = 0;         % Number of parallel workers used

save_to_file = false;
results_file = sprintf('D%i%i_S%i_Random_Walk', n,k,d);
image_file = results_file;
temp_file = 'parsave';

% -------------------------------------------------------------------------
% Set up the parallel pool
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

% -------------------------------------------------------------------------
% Set up variables for both the parallel and sequential parts
% -------------------------------------------------------------------------
% Parallel variables
confs_0 = zeros(n, d+1, nReps);
dms_0 = zeros(n,n, nReps);
bd_times_0 = zeros(nReps,2);

% Create a parallel.pool.Constant variable. The workers will use them to
% save results while the program runs. In case the program crashes, these
% files can be used to read intermediate results
spmd
    t = getCurrentTask();
    myFname = sprintf('temp/checkpoint_%s_%02d', temp_file, t.ID); % each worker gets a unique filename
    myMatfile = matfile(myFname, 'Writable', true);
    
    % Seed the variables in the matfile object
    myMatfile.confs_0 = confs_0;
    myMatfile.dms_0 = dms_0;
    myMatfile.bd_times_0 = bd_times_0;
    myMatfile.gotResult = false(1,nReps);
end
myMatfileConstant = parallel.pool.Constant(myMatfile);

% Sequential variables
myFname = 'temp/checkpoint_non_parallel';
myMatfile_seq = matfile(myFname, 'Writable', true);

confs_n = zeros(n,d+1,nSteps);
dms_n = zeros(n,n,nSteps);
bd_times_n = zeros(nSteps,2);

myMatfile_seq.confs_n = confs_n;
myMatfile_seq.dms_n = dms_n;
myMatfile_seq.bd_times_n = bd_times_n;

t_setup = toc;

% -------------------------------------------------------------------------
% Sample uniformly at random
% -------------------------------------------------------------------------
% Record when the simulation itself started
tic

parfor i=1:nReps
    % Sample n points from S^d
    X = randn(n,d+1);
    norm = sqrt(sum(X.^2,2));
    X = X./repmat(norm, [1,d+1]);
    
    % Get the distance matrix
    dm = L2_distance(X',X');
    dm = 2*asin(dm/2);

    % Compute the Vietoris-Rips persistence diagrams of the metric space
    % with distance matrix dm
    [tb,td] = persistence_matrix(dm);

    % Save results
    confs_0(:,:,i) = X;
    dms_0(:,:,i) = dm;
    bd_times_0(i,:) = [tb,td];
    
    % And save them to disk (for redundancy)
    matfileObj = myMatfileConstant.Value;
    matfileObj.confs_0(:,:,i) = X;
    matfileObj.dms_0(:,:,i) = dm;
    matfileObj.bd_times_0(i,:) = [tb,td];
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
t_par = toc;

delete(gcp('nocreate'));

% -------------------------------------------------------------------------
% Execute the random walk
% -------------------------------------------------------------------------
% Choose a configuration to start from:
% The cross-polytope
X0 = eye(d+1);
X0 = X0(1:n/2,:);
X0 = [X0; -X0];

tic
i = 1; bad_j=0; backtrack = 0;
while i <= nSteps
    % Perturb the configuration after the first iteration
    if i==1
        X=X0;
    else
        X = perturb_config_Sn(X, sigma);
    end

    % Get the distance matrix
    dm = L2_distance(X',X');
    dm = 2*asin(dm/2);

    % Compute persistence
    [tb,td] = persistence_matrix(dm);

    if tb+td ~= 0
        % Now we check if it is close to known persistence
        % Calculate distance to known configurations
        post_s = [bd_times_0; bd_times_n(1:(i-1-backtrack),:)];
        nbhd_s = d_bottleneck([tb,td], post_s);
        nbhd_density_s = sum(nbhd_s < eps)/size(post_s,1);
        
        post_pre = post_s(1:(end-1),:);
        nbhd_pre = d_bottleneck(post_s(end,:), post_pre);
        nbhd_density_pre = sum(nbhd_pre < eps)/size(post_pre,1);
        
        %alpha = nbhd_density_s/nbhd_density_pre;    % try going to higher density
        % try going to lower density
        if nbhd_density_s ~= 0
            alpha = nbhd_density_pre/nbhd_density_s;
        else
            % If nbhd_density==0, there are no points around [tb,td] so we
            % definitely want to go there
            alpha = 1;
        end
        
        % Display density
        fprintf('%0.3d, ', alpha);
        
        % Decide to accept by tossing a coin
        if rand <= alpha
            % Add them to the list
            confs_n(:,:,i) = X;
            dms_n(:,:,i) = dm;
            bd_times_n(i,:) = [tb,td];
            
            % Save results to disk every once in a while
            % (around every 10% of the repetitions)
            m = ceil(nSteps/10);
            if rem(i, m)==0
                i0 = i-m+1;
                myMatfile_seq.confs_n(:,:,i0:i) = confs_n(:,:,i0:i);
                myMatfile_seq.dms_n(:,:,i0:i) = dms_n(:,:,i0:i);
                myMatfile_seq.bd_times_n(i0:i,:) = bd_times_n(i0:i,:);
            end
            
            fprintf('\n%1d\n-------------\n',i);
            disp(i);
            i = i +1;   % advance the counter
            bad_j = 0;
            backtrack = 0;
        else
            % If we don't accept, we perturb the previous configuration
            if i==1
                X = X0;
            else
                X = confs_n(:,:,i-1);
            end
            % We will stop searching this configuration if it produces too
            % many useless points
            fprintf("(bad %d), ", bad_j);
            bad_j = bad_j + 1;
            if rem(bad_j,100)==0
                fprintf('\n');
            end
            if bad_j == 50
                bad_j = 0;
                backtrack = min(i-1, backtrack+1);
                X = confs_n(:,:,i-backtrack);
                fprintf('\n');
            end
        end
    end
    % If there is no persistence, bd_times remains 0 at that row
end % Of one of the nReps repetitions

t_seq = toc;
fprintf('Preparations: %03f\n',t_setup);
fprintf('Parallel calculations: %0.3f\n',t_par);
fprintf('Sequential calculations: %0.3f\n',t_seq);

%% Plot the results
figure();
hold on;
plot(bd_times_n(:,1), bd_times_n(:,2), 'r.')
plot(bd_times_0(:,1), bd_times_0(:,2), 'b.')
plot_reference_lines();
title(sprintf('D_{%i,%i}(S^%i)', n,k,d));
hold off;

%% Save the results
if save_to_file
    % Put both lists together
    confs = cat(3,confs_0,confs_n); %#ok<UNRCH>
    dms = cat(3,dms_0,dms_n);
    bd_times = cat(1,bd_times_0,bd_times_n);

    exclude = 't|p|myFname|myMatfile|myMatfileConstant|dX_constant|myMatfile_seq|exclude';
    save(sprintf('results/%s.mat', results_file), '-regexp', sprintf('^(?!(%s)$).', exclude));
    % save(sprintf('results/%s.mat', results_file));
    saveas(gcf, sprintf('results/%s.png', image_file));
end
