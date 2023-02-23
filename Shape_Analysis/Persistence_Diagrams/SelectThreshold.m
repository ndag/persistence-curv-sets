% Given a distance matrix dm, this function finds the entry d of dm such
% that the percentage of distances strictly less than d is smaller than s
% Notice 0 < s <= 1
function t = SelectThreshold(dm, s)
% Use squareform if we have the square matrix
if size(dm,1) > 1
    dm = squareform(dm);
end

% Number of possible edges
N = length(dm);

% Sort the distances
d_sort = sort(dm);

% Find the index idx such that i/N < s
idx = floor(N*s);

% Verify if t=d_sort(idx) satisfies what we want
% If it doesn't, then t is the largest distance smaller than d_sort(idx)
t = d_sort(idx);
if sum(d_sort<=t)/N > s
    if sum(d_sort<t) > 0
        t = max( d_sort(d_sort<t) );
    else
        disp('Warning. SelectThreshold returned 0')
        t=0;
    end
end

end