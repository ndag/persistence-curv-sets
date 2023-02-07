function X = regular_polygon(numSides, theta0)
    if nargin<2
        theta0 = pi/2;
    end
    angles = (0:numSides)*(2*pi/numSides);
    x = cos(angles+theta0);
    y = sin(angles+theta0);
    
    X = [x',y'];
end
