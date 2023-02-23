% Calculates the bottleneck distance between the one-point diagrams a in A
% and b in B, for every a, b in A.
% Each row of A and B should be a persistence diagram. In other words:
% - A = Nx2 matrix
% - B = Mx2 matrix
% - dB= NxM matrix where dB(i,j) is the distance between A(i,:) and B(j,:)
function dB = d_bottleneck(A,B)
    % Calculate the L_infinity distance
    difference = abs(permute(A, [1,3,2]) - permute(B, [3,1,2]));
    L_inf = max(difference,[], 3);
    
    % Calculate the persistence of all points
    pers_A = A(:,2)-A(:,1);
    pers_B = B(:,2)-B(:,1);
    
    % Find the maximum persistence of each pair
    pers = max(pers_A,pers_B');
    
    % Bottleneck distance between one-point diagrams is the minimum of
    % either the L_infinity distance or half the persistence.
    dB = min(L_inf, 0.5*pers);
end