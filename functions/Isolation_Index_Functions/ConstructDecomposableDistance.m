% Construct a symmetric matrix with zeros on the diagonal as described on
% the paper "A Canonical Decomposition Theory for Metrics on a Finite Set".
% In some cases, the result is actually a distance matrix.

% INPUT:
% % alpha: A vector with n positive components. These are the isolation
% indices.
% % splits: A cell array of size [n,2] of row vectors with positive
% integers such that the union of every row splits{i,1} U splits{i,2}]
% equals the vector 1:N for some N>0. These are the dm-splits.

% alpha(i) is the isolation index of the split splits{i,1}, splits{i,2}.
function dm = ConstructDecomposableDistance(alphas, splits)
    % Number of points in the resulting metric space
    N = numel(splits{1,1})+numel(splits{1,2});
    
    % Number of summands
    n = size(splits,1);
    
    % Add all the split-metric summands
    dm = zeros(N,N);
    for i=1:n
        % Construct the split metric d_A,B
        A = splits{i,1};
        B = splits{i,2};
        dm_s = ConstructSplitMetric(A,B);
        
        % Add the summand
        dm = dm+alphas(i)*dm_s;
    end
end