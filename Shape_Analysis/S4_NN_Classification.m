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
classes = cell(size(names));
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

% Load the names of the distance matrices and repeat the experiment for
% each 
matrix_names = dir('Results/Distances/*.mat');
matrix_names = {matrix_names.name};
matrix_names = cellfun(@(s) strsplit(s, '.'), matrix_names, 'UniformOutput', false);
matrix_names = cellfun(@(s) s{1}, matrix_names, 'UniformOutput', false);

for idm=1:numel(matrix_names)
    matrix_name = matrix_names{idm};

    dm = load(sprintf('Results/Distances/%s.mat', matrix_name));
    if strcmp(matrix_name(1:2), 'dH')
        dm = dm.dH;
    elseif strcmp(matrix_name(1:2), 'dW')
        dm = dm.dW;
    elseif strcmp(matrix_name(1:2), 'dB')
        dm = dm.dB;
    end
    dm = squareform(dm);
    dm = dm(1:nShapes, 1:nShapes);

    cm_all = zeros(nLabels, nLabels, nTests);
    Reps = zeros(nTests, nLabels);
    for idt=1:nTests
        % Choose one representative for each label
        for idl=1:nLabels
            candidates = find(classes==idl);
            Reps(idt, idl) = candidates( randi(length(candidates)) );
        end

        % Find the nearest neighbor for all points
        dmt = dm(Reps(idt,:), :);
        [~, classes_predicted] = min(dmt);      % I contains the row index of the minimum in each column

        % Syntax: confusionMat(known, predicted)
        % Output: Row i contains the predicted labels of all entries with true
        % class i
        cm_all(:,:,idt) = confusionmat(classes, classes_predicted);
    end

    % Results
    % Compile experiments
    cm = sum(cm_all,3);
    f = confusionchart(cm, labels, 'Normalization', 'row-normalized');

    % Compute percentage of missclassification
    miss = sum( cm-diag(diag(cm)), 'all')/(nTests*nShapes);

    f.Title = sprintf('%s (error: %0.2f%%)',...
        strrep(matrix_name, '_', '\_'),...
        miss*100);

    export_fig(sprintf('Results/Distances/%s_confusion.pdf', matrix_name));
end
