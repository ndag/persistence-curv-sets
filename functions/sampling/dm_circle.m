% Returns the distance matrix of a subset of the circle with geodesic
% distance. We model the circle as a quotient [0,2*R] / (0 ~ 2R), and
% geodesic distance
% d(x,y) = min(|x-y|, 2R-|x-y|)
% - X: Vector of points in [0,2R]
% - R: Diameter of the circle (so that the perimeter is 2*R)
function dm = dm_circle(X,R)
    if nargin<2
        R=pi;
    end
    
    diff = abs(X-X');
    
    dm = min(diff, 2*R-diff);