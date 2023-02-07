% The function simple_polygon constructs the vertices [x,y] of a random
% simple polygon. Here, I put approximately N points on the sides of the
% polygon as uniformly distributed as I can.
% For that, I compute the length of all sides and I insert an
% equidistributed amount of points in each side proportional to its length.
function [xx, yy] = simple_curve(X,Y, N)
    % Number of sides in the polygon
    n = numel(X);
    
    % Compute the length of the sides
    lengths = zeros(n-1,1);
    for i=1:(n-1)
        lengths(i) = sqrt((X(i+1)-X(i))^2 + (Y(i+1)-Y(i))^2);
    end
    % Note: the last point always equals the first
    
    % Total length
    L = sum(lengths);
    
    % The number of points in each side is proportional to the length of
    % the side
    Ns = round(N*lengths/L);
    
    % Create the intermediate points in each side
    xx = [];
    yy = [];
    for i=1:(n-1)
        x_row = linspace(X(i),X(i+1),Ns(i));
        y_row = linspace(Y(i),Y(i+1),Ns(i));
        
        xx = [xx, x_row(1:end-1)]; %#ok<*AGROW>
        yy = [yy, y_row(1:end-1)];
    end
end


% % Number of sides in the polygon
% m = numel(X);
% 
% % Compute the length of the sides
% lengths_2 = zeros(m-1,1);
% for i=1:(m-1)
%     lengths_2(i) = sqrt((X(i+1)-X(i))^2 + (Y(i+1)-Y(i))^2);
% end