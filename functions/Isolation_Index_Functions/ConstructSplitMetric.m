% This function constructs the split-metric d_A,B given a split [A,B].
% Here,
% d_A,B(x,y)=1 iff x \in A, y \in B (or viceversa) and
% d_A,B(x,y)=0 iff x,y \in A or x,y \in B.
% The input are 2 row vectors A,B with integer entries. We want both A and
% B to have distinct entries, both within themselves, as well as between
% each other.
function dm = ConstructSplitMetric(A, B)
    n_a = numel(A);
    n_b = numel(B);
    
    dm = zeros(n_a+n_b);
    dm(A,B) = 1;
    dm(B,A) = 1;
end