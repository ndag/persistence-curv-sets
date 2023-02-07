function [] = plot_reference_lines(r)
    if nargin == 0
        r = pi;
    end
    
    % D_{4,1}(S1)
    xx1 = 0:0.001:r;
    xx2 = 0.5*r:0.001:r;
    
    plot(xx1,xx1,'red');                % diagonal
    plot(xx2,xx2,'red');
    plot(xx2,ones(size(xx2))*r,'red');    % top line

    % Left "diagonal"
    tb = sqrt(2):0.001:sqrt(8/3);
    td = sqrt(4*tb.^2.*(3-tb.^2)./(4-tb.^2));
    
    % The description above is in Euclidean distance. We need to convert it
    % to Spherical distance
    plot(2*asin(tb/2), 2*asin(td/2), 'red');
    
    % Persistence of the regular hexagon slided to the north pole
    tt = 0:0.001:1;
    xx = acos(1-1.5*tt.^2)*r/pi;
    yy = acos(1-2*tt.^2)*r/pi;
    plot(xx,yy,'red');

    % Left boundary of D_6,2(S^1)
    tt = 2/3*r:0.01:r;
    plot(2/3*r*ones(size(tt)), tt, 'red')
    
    % Persistence of the square with two extra points slided to the north pole
    % tt = 0:0.001:1;
    % xx = acos(1-(1-cos(3*pi/4))*tt.^2)*r/pi;
    % yy = acos(1-2*tt.^2)*r/pi;
    % plot(xx,yy,'red');
end