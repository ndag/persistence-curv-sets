clear; close all; clc;

tic
save_to_file = false;
results_file = 'D62_S2_Random_Walk';
image_file = 'D62_S2_Random_Walk';

% Load functions and previous configurations
addpath('../../functions/persistence/');
addpath('../../functions/functions_S2/');
load('results/2021_06_22/D62_S2.mat')

% Change names
confs_0 = confs;
dms_0 = dms;
bd_times_0 = bd_times;

% Parameters for the simulation
n = 6;              % Number of points in each configuration
d = 2;              % Dimension of the sphere
sigma = 0.075;      % Standard deviation of the normal perturbation
eps = 0.005;        % Radius of density neighborhood
density = 0.01;     % Maximum density
rng(1);             % set random seed

% Number of configurations to sample
nSteps = 1000;

% List to save the results
confs_n = zeros(n,d+1,nSteps);
dms_n = zeros(n,n,nSteps);
bd_times_n = zeros(nSteps,2);

% Initial configuration
% Regular polytope
X0 = [1,0,0; -1,0,0; 0,+1,0; 0,-1,0; 0,0,+1; 0,0,-1];

% -------------------------------------------------------------------------
% Repeat the experiment
% -------------------------------------------------------------------------
% Plot all current points
figure();
hold on;
plot(bd_times_0(:,1), bd_times_0(:,2), 'b.')
plot_reference_lines();
drawnow;

t1 = toc;
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

    confs_n(:,:,i) = X;
    dms_n(:,:,i) = dm;

    % Create a metric space with the distance matrix of the 6 points in L
    % and compute its Vietoris-Rips complex
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
            bd_times_n(i,:) = [tb,td];
            confs_n(:,:,i) = X;
            dms_n(:,:,i) = dm;
            
            % Plot the new point in green, and paint the previous point in
            % red
            plot(bd_times_n(i,1), bd_times_n(i,2), 'g.');
            if i>1
                plot(bd_times_n(i-1,1), bd_times_n(i-1,2), 'r.');
            end
            drawnow;
            % pause;
            
            fprintf("\n%1d\n-------------\n",i);
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
                fprintf("\n");
            end
            if bad_j == 50
                bad_j = 0;
                backtrack = min(i-1, backtrack+1);
                X = confs_n(:,:,i-backtrack);
                fprintf("\n");
            end
        end
    end
    % If there is no persistence, bd_times remains 0 at that row
end % Of one of the nReps repetitions

t = toc;
fprintf("Preparations: %03f\n",t1);
fprintf("Calculations: %0.3f\n",t);

% Plot the results
figure();
hold on;
plot(bd_times_0(:,1), bd_times_0(:,2), 'b.')
plot(bd_times_n(:,1), bd_times_n(:,2), 'r.')
plot_reference_lines();
hold off;

% Save the results
if save_to_file
    % Put both lists together
    confs = cat(3,confs_0,confs_n); %#ok<UNRCH>
    dms = cat(3,dms_0,dms_n);
    bd_times = cat(1,bd_times_0,bd_times_n);
    
    % exclude = 'myFname|myMatfile|exclude';
    % save(sprintf('./results/%s.mat', results_file), '-regexp', sprintf('^(?!(%s)$).', exclude));
    save(sprintf('./results/%s.mat', results_file));
    saveas(gcf, sprintf('%s.png', image_file));
end
