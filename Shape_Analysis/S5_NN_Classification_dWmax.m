clear; close all; clc; 

% Use 1-NN to classify the shapes
nTests = 2000;
rng(1)

% Workaround for a weird error when using confusionchart
f = figure('visible', 'off');
set(f, 'Color', 'w')
hold off

% Load the name of the shapes and determine the class they belong to
names = dir('models/processed/*.mat');
names = {names.name};
nShapes = numel(names);

% The class labels are obtained by keeping only the first part of the
% string (before any dashes, dots or numbers)
classes = cell(1,nShapes);
for idx=1:nShapes
    name = names{idx};
    letters = isletter(name);

    % Get the index of the first non-letter and keep everything before it
    i0 = find(letters==0);
    i0 = i0(1);

    classes{idx} = name(1:i0-1);
end

labels = unique(classes);
nLabels = numel(labels);

% Turn classes into numbers
classes = cellfun(@(s) find(cellfun(@(cls) strcmp(s, cls), labels)), classes);

% ----------------------
% Classification labels
% ----------------------
% Compute middle point of each class (to label the distance matrix heatmap)
C1 = arrayfun( @(n) find(classes==n, 1, 'first'), 1:nLabels);
C2 = arrayfun( @(n) find(classes==n, 1,  'last'), 1:nLabels);
labCenter = round((C1+C2)/2);

% Load the names of the distance matrices
matrix_names = dir('Results/Distances/dW*.mat');
matrix_names = {matrix_names.name};
matrix_names = cellfun(@(s) strsplit(s, '.'), matrix_names, 'UniformOutput', false);
matrix_names = cellfun(@(s) s{1}, matrix_names, 'UniformOutput', false);

% Load distance matrices
dms = zeros(nShapes,nShapes, numel(matrix_names));
for idm=1:numel(matrix_names)
    matrix_name = matrix_names{idm};

    file = load(sprintf('Results/Distances/%s.mat', matrix_name));
    if strcmp(matrix_name(1:2), 'dH')
        dm = file.dH;
    elseif strcmp(matrix_name(1:2), 'dW')
        dm = file.dW;
    elseif strcmp(matrix_name(1:2), 'dB')
        dm = file.dB;
    end
    % Turn Inf into NaN
    dm(isinf(dm)) = nan;
    dm = squareform(dm);
    dms(:,:,idm) = dm(1:nShapes,1:nShapes);
end

% Function to reweight the distance matrices
wmax = @(M, w) max(M.*permute(w, [1,3,2]), [], 3);

% Use the same representatives in all trials
Reps = zeros(nTests, nLabels);
for idt=1:nTests
    for idl=1:nLabels
        candidates = find(classes==idl);
        Reps(idt, idl) = candidates( randi(length(candidates)) );
    end
end

% Function to get the error from NN classification
NN = @(w) NN1_function( wmax(dms,w), classes, Reps);

% Optimize NN
w0 = 1/nanmax(dms, [], [1,2]);
w0 = permute(w0, [1,3,2]);
% w0(2) = 0; % Modify w0
[w, fval] = fminsearch(NN, w0, optimset('Display', 'iter'));

% Re-do the experiment to obtain the confusion matrix
disp(w0)
disp(w)
[miss, cm] = NN(w);

%%
% ------------------------
% Plot and save results
% ------------------------
% Create a new file name each time we run this script (to not overwrite
% previous results)
idx=0;
file_exists = true;
while file_exists
    % If the file already exists, we check the next number
    idx = idx+1;
    file_name = sprintf('Results/Distances/optim/dW_(%0.1f)_%i', 100*miss, idx);
    file_exists = isfile(sprintf('%s.png', file_name));
end % And repeat until we get a new file

% Weighted distance matrix
f = figure;

colormap('hot');
imagesc(wmax(dms,w));
hold on
title(sprintf('dW (error:%0.2f%%)\n %0.1f, %0.1f, %0.1f', 100*miss, w(1), w(2), w(3)))
colorbar;

axis('equal')
xlim([1,nShapes]);
ylim([1,nShapes]);

set(gca, 'xtick', labCenter, 'xticklabel', labels)
set(gca, 'ytick', labCenter, 'yticklabel', labels)

% Stuff to remove whitespace
set(f, 'Color', 'w')
export_fig(sprintf('%s.pdf', file_name));
hold off;

% Confusion matrix
f = confusionchart(cm, labels, 'Normalization', 'row-normalized');

f.Title = sprintf('dW (error:%0.2f%%)\n %0.1f, %0.1f, %0.1f',...
                   100*miss, w(1), w(2), w(3));
export_fig(sprintf('%s_Conf.pdf', file_name));

% Save results to file
% mFile = matfile(sprintf('%s.mat', file_name));
% mFile.dms = dms;
% mFile.w0 = w0;
% mFile.w = w;
% mFile.cm = cm;
% mFile.miss = miss;
