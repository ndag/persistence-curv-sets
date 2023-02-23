clear; close all; clc;

% Load packages
addpath('../functions/persistence')
addpath('../functions/distances')
addpath('../functions/sampling')
addpath('../functions/export_fig')
addpath('../functions/bottleneck')

addpath('Persistence_Diagrams')

% Parameters
maxHomDim = 2;  % Maximum Homology dimension to compute with Ripser

% Choose the largest nL so that Ripser can compute PH_k with k=maxHomDim
nL = 500;       % Number of landmark points

% Percentege of edges that we want to allow (see function SelectThreshold)
s_edges = 0.07;

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

%% 
% Subsample
P1_SubsampleMMSpaces
disp('------------------------------ ')
disp('Subsampling done')
disp('------------------------------ ')

% Compute persistence diagrams
S1_ComputePH
disp('------------------------------ ')
disp('PH Done')
disp('------------------------------ ')

S2_Compute_dB
disp('------------------------------ ')
disp('dB Done')
disp('------------------------------ ')

[labels, labCenter] = LabelsForGraphs(names);
S3_Show_dB
disp('------------------------------ ')
disp('Plot dB Done')
disp('------------------------------ ')

% Clean temp files
which_dir = 'temp';
dinfo = dir(which_dir);
dinfo([dinfo.isdir]) = [];   %skip directories
filenames = fullfile(which_dir, {dinfo.name});
delete( filenames{:} )
