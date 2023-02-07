function [alpha, A] = IsolationIndexMatrixOrder2(dm)
    n = size(dm,1);
    
    alpha = zeros(n,1);
    A = zeros(n,n);
    for i=1:n
        alpha(i) = IsolationIndex(dm, i);
        % A(i,i) = alpha(i);
        
        for j=i+1:n
            index = IsolationIndex(dm, [i,j]);
            A(i,j) = index;
            A(j,i) = index;
        end
    end
end