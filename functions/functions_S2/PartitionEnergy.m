% Given a configuration X in S^2, we calculate the functional K(X) in
% equation (3.4) of Centroidal Voronoi Tessellations: Applications and
% Algorithms. Recall that V is the Voronoi partition of the sphere
% generated by X.
function K = PartitionEnergy(X)
    n = size(X,1);
    Y = permute(X,[3,4,1,2]);
    
    % The Jacobian in spherical coordinates
    J = @(s,t) sin(t);
    
    % A point in S^2 given in spherical coordinates
    P = @(s,t) repmat(...
        cat(4, cos(s).*sin(t), sin(s).*sin(t), cos(t)),...
        1,1,n);
    
    % The potiential of a point is the squared distance to the closest
    % point.
    F = @(s,t) min(sum((Y-P(s,t)).^2,4), [], 3);
    
    % The paper split the integral into the Voronoi cells. However, the
    % Voronoi cells contain the points closest to its center than any other
    % center, so we can instead integrate F as above over the whole sphere
    % (multiplied by the Jacobian).
    K = integral2( @(x,y) F(x,y).*J(x,y), 0,2*pi,0,pi);
end