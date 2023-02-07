% Constructs a regular n-simplex inscribed on the sphere S^(n-1).
function [X, dm] = build_regular_simplex(n)
    if n==1
        % When we reach n==1, S^0 is just two points at distance pi
        X = [1; -1];
        
        dm = [0, 2; 2, 0];
        dm = 2*asin(dm/2);
    else
        % The first point is the north pole
        X_1 = [zeros(1,n-1), 1];

        % We construct a regular simplex in S^(n-1)
        X_next = build_regular_simplex(n-1);
        X_next = [sqrt(1-1/n^2)*X_next,...
                 -1/n*ones(size(X_next,1),1)];

        % And rescale to fit it in S^n
        X = [X_1; X_next];

        dm = L2_distance(X',X');
        dm = 2*asin(dm/2);
    end
end
