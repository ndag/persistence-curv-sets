function [tb, td] = bd_times_matrix(M)
    % Takes an n-by-n matrix M with n=2k+2 and calculates t_b^(M) and
    % t_d^(M) as given by my algorithm.
    % This is NOT the same as persistent homology because tb can be greater
    % than td.
    t_bd = maxk(M,2);
    td = min(t_bd(1,:));
    tb = max(t_bd(2,:));
end
