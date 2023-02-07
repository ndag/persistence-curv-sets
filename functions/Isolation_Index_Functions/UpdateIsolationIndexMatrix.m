% This function receives a distance matrix dm, a list of partial dm-splits
% (splits), a vector with their isolation indices (alpha), and a new point
% x_new. All partial splits A,B are partitions of a fixed set Y = A \cup B.
% We will calculate the isolation index in Y \cup {x_new} of:
% - {x_new}, Y
% - A U {x_new}, B
% - A, B U {x_new},
% where A,B run over all splits in the array splits. We only keep the new
% splits that have non-zero isolation index, and return the results.
function [alpha_new, splits_new] = UpdateIsolationIndexMatrix(x_new, splits, alpha, dm)
    Y = [splits{1,1}, splits{1,2}];
    
    N_splits = size(splits,1);
    alpha_new = zeros(N_splits*2+1, 1);
    splits_new = cell(N_splits*2+1, 1);
    
    % Check if {x_new}, Y has non-zero isolation index
    n_Y = numel(Y);
    Y_new = [Y, x_new];
    dm_new = dm(Y_new, Y_new);
    alpha_Y = IsolationIndex(dm_new, n_Y+1, false);
    % If it does, we add it to the new list
    if alpha_Y > 0
        alpha_new(1) = alpha_Y;
        splits_new{1,1} = Y;
        splits_new{1,2} = x_new;
    end
    
    % Do the same for all other non-trivial splits
    for i=1:N_splits
        alpha_0 = alpha(i);
        A = splits{i,1};
        B = splits{i,2};
        
        % See if A U {x_new}, B and A, B U {x_new} have non-zero isolation
        % index
        alpha_A = UpdateIsolationIndex(x_new, A, B, alpha_0, dm);
        alpha_B = UpdateIsolationIndex(x_new, B, A, alpha_0, dm);
        
        % If they do, we add them to the new list
        if alpha_A > 0
            alpha_new(2*i) = alpha_A;
            splits_new{2*i,1} = [A, x_new];
            splits_new{2*i,2} = B;
        end
        if alpha_B > 0
            alpha_new(2*i+1) = alpha_B;
            splits_new{2*i+1,1} = A;
            splits_new{2*i+1,2} = [B, x_new];
        end
    end
    
    % Finally, we filter out the splits whose isolation index is 0
    non_splits = find(alpha_new == 0);
    alpha_new(non_splits) = [];
    splits_new(non_splits,:) = [];
end