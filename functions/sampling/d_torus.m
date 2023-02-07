function d = d_torus(X, Y, p, R1, R2)
    % Modeling the torus as the square [0,R1]*[0,R2] with the usual
    % identifications, equiv. (S^1)x(S^1), the metric here is the Lp
    % product of the metrics of S^1. Equivalently, under this
    % parametrization, the distance in the torus is the usual Lp distance
    % of R^2.
    % X and Y are n-by-2 matrices.
    if nargin < 3
        p=2;
    end
    if nargin < 5
        R1=pi;
        R2=pi;
    end
    % Use the metric on the circle
    d1 = X(:,1)-Y(:,1)';
    d2 = X(:,2)-Y(:,2)';
    
    d1 = min(abs(d1), 2*R1-abs(d1));
    d2 = min(abs(d2), 2*R2-abs(d2));
    
    %L_infty product
    if isinf(p)
        d = max(d1,d2);
    % Lp product
    else
        d = nthroot(d1.^p+d2.^p, p.*ones(size(X,1)) );
    end
end