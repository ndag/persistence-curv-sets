function [miss, cm] = NN1_function(dm, classes, Reps)
nShapes = size(dm,1);

nTests = size(Reps,1);
nLabels = size(Reps,2);

cm_all = zeros(nLabels, nLabels, nTests);
for idt=1:nTests
    % Find the nearest neighbor for all points
    dmt = dm(Reps(idt,:), :);
    [~, classes_predicted] = min(dmt);      % I contains the row index of the minimum in each column

    % Syntax: confusionMat(known, predicted)
    % Output: Row i contains the predicted labels of all entries with true
    % class i
    cm_all(:,:,idt) = confusionmat(classes, classes_predicted);
end

% Compile results of all tests
cm = sum(cm_all,3);

% Compute percentage of missclassification
miss = sum( cm-diag(diag(cm)), 'all')/(nTests*nShapes);
