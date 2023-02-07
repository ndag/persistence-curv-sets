function [tb,td] = persistence_matrix(M)
    % Calculates the persistent homology of a metric space with distance
    % matrix M using bd_times_matrix. Since that function can return tb>td,
    % here we check for that condition and return td=tb=0 if that's true.
    [tb,td] = bd_times_matrix(M);
    
    % If td<=tb, then there is no live cycle.
    if td<=tb
        tb=0;
        td=0;
    end
    % But if tb<td, we found the persistent homology so we return those
    % values without change.
    