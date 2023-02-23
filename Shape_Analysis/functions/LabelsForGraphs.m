% ----------------------
% Classification labels
% ----------------------
function [labels, labCenter] = LabelsForGraphs(names)
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

% Compute middle point of each class (to label the distance matrix heatmap)
C1 = arrayfun( @(n) find(classes==n, 1, 'first'), 1:nLabels);
C2 = arrayfun( @(n) find(classes==n, 1,  'last'), 1:nLabels);
labCenter = round((C1+C2)/2);

end