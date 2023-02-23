clear; close all; clc;

% Load packages
addpath('functions')
addpath('../functions/distances')
addpath('../functions/sampling')
addpath('../functions/readObj/')

% Create folders to save results
if ~exist('temp', 'dir')
    mkdir('temp')
end

input_folder = 'models/processed';
if ~exist(input_folder, 'dir')
    mkdir(input_folder)
end

param_folder = 'Results/Parameters';
if ~exist(param_folder, 'dir')
    mkdir(param_folder)
end

% List of names of shapes
names = dir('models/*.mat');
names = {names.name};
names = cellfun(@(s) strsplit(s, '.'), names, 'UniformOutput', false);
names = cellfun(@(s) s{1}, names, 'UniformOutput', false);

% Compute histogram centers
nBins0 = 850;       % Number of histogram centers in dimension 0
nBins1 = 946;       % Approximate number of centers in dimensions > 0

P0_Histograms_Centers
disp(' ------------------------------ ')
disp('Histogram centers Done')
disp(' ------------------------------ ')

% Compute distance matrices
% NFPS = 4000;        % Number of points to subsample with FPS
NFPS = 1000;        % Number of points to subsample with FPS

rng(1);
P1_Compute_dm_Sumner
disp(' ------------------------------ ')
disp('DMs Done')
disp(' ------------------------------ ')
