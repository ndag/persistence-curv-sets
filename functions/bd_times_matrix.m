function [tb, td] = bd_times_matrix(M)
    % Takes an n-by-n matrix M with n=2k+2 and calculates t_b^(M) and
    % t_d^(M) as given by my algorithm.
    % This is NOT the same as persistent homology because tb can be greater
    % than td.
    n = size(M,1);
    
    for i=1:n
        M(i,:) = sort(M(i,:));
    end
    
    tb = max(M(:,n-1));
    td = min(M(:,n));
end