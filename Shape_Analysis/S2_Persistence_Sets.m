clear; close all; clc;

% Load packages
addpath('../functions/persistence')
addpath('../functions/distances')
addpath('../functions/sampling')
addpath('../functions/export_fig')
addpath('../functions/emd')

addpath('Persistence_Sets')

% Parameters
Ks = [0, 1,2];              % Dimensions to calculate principal persistence sets
nRepsK = [0, 10^4, 10^4];   % Number of samples to take in each persistence set
% Note: Leave a placeholder for dimension 0 in nRepsK. We don't need to cap
% nReps in dimension 0, but we still need a value there

% Create folders to save results
input_folder = 'models/processed';

results_folder = 'Results';
path = results_folder;
if ~exist(path, 'dir')
    mkdir(path)
end

path = sprintf('%s/Distances', results_folder);
if ~exist(path, 'dir')
    mkdir(path)
end

% List of names of shapes
names = dir(sprintf('%s/*.mat', input_folder));
names = {names.name};
names = cellfun(@(s) strsplit(s, '.'), names, 'UniformOutput', false);
names = cellfun(@(s) s{1}, names, 'UniformOutput', false);

useShortlist = false;
if useShortlist
    % load('Results/Parameters/Shortlist.mat');
    short = 1:3;
    names = names(short);
end
nShapes = numel(names);

%% Repeat the experiment for each k
for ik=1:numel(Ks)
    k = Ks(ik);
    n = 2*k+2;
    nReps = nRepsK(ik);
    
    % Create folder to save results
    path = sprintf('%s/D_%i_%i', results_folder,n,k);
    if ~exist(path, 'dir')
        mkdir(path)
    end
    
    % Load histogram parameters
    if k==0
        Params = matfile('Results/Parameters/Hist_Dn0.mat');
    else
        Params = matfile('Results/Parameters/Hist_947.mat');
    end
    
    disp('------------------------------ ')
    disp('------------------------------ ')
    fprintf(' START k=%i\n', k)
    disp('------------------------------ ')
    disp('------------------------------ ')
    
    % Call the scripts
    rng(1);
    if k==0
        S1_Compute_PS0
    else
        S1_Compute_PS
    end
    disp('------------------------------ ')
    disp('PS Done')
    disp('------------------------------ ')

    rng(1);
    S2_Compute_Histogram_Vor
    disp('------------------------------ ')
    disp('Histogram Done')
    disp('------------------------------ ')

    rng(1);
    if k==0
        S3_dH_PS0
    else
        S3_dH_PS
    end
    disp('------------------------------ ')
    disp('dHausdorff done')
    disp('------------------------------ ')
    
    rng(1);
    S4_dW_Histogram_Vor
    disp('------------------------------ ')
    disp('dWasserstein done')
    disp('------------------------------ ')
end

[labels, labCenter] = LabelsForGraphs(names);
SE1_Show_dH
disp('------------------------------ ')
disp('Plot dH done')
disp('------------------------------ ')

SE2_Show_dW
disp('------------------------------ ')
disp('Plot dW done')
disp('------------------------------ ')

SE3_Max_Dists
disp('------------------------------ ')
disp('Plot d-Max done')
disp('------------------------------ ')

% Clean temp files
which_dir = 'temp';
dinfo = dir(which_dir);
dinfo([dinfo.isdir]) = [];   %skip directories
filenames = fullfile(which_dir, {dinfo.name});
delete( filenames{:} )
