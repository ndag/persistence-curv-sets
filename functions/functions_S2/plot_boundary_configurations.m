function [] = plot_boundary_configurations(color)
if nargin<1
    color='red';
end

numH = 25;
start = pi/5;
HH = start:((pi/2-start)/(numH-1)):pi/2;
T = pi/3;

bd_times = zeros(numH,2);

for i=1:numH
    H = HH(i);
    theta1 = pi/2-H;

    % Solve a polynomial to find H2 such that the second inter-triangle
    % distance equals the side of the larger triangle
    A1 = cos(theta1);
    A2 = cos(T+4*pi/3)*sin(theta1);
    R = roots([9, -3*A1, -5, -A1]);
    R = R(IsNear(imag(R), 0, 10^(-5)));      % Filter out complex roots
    H2 = max(acos(real(R)))-pi/2;

    % Equilateral triangle
    t1 = 0;
    t2 = 2*pi/3;
    t3 = 4*pi/3;

    % A shifted equilateral triangle
    t4 = t1+T;
    t5 = t2+T;
    t6 = t3+T;

    % The first triangle is at height H above the equator
    s1 = pi/2-H;
    s2 = pi/2-H;
    s3 = pi/2-H;

    % The second triangle is below the first one
    s4 = pi/2+H2;
    s5 = pi/2+H2;
    s6 = pi/2+H2;

    X = [cos(t1)*sin(s1), sin(t1)*sin(s1), cos(s1);
       cos(t2)*sin(s2), sin(t2)*sin(s2), cos(s2);
       cos(t3)*sin(s3), sin(t3)*sin(s3), cos(s3);
       cos(t4)*sin(s4), sin(t4)*sin(s4), cos(s4);
       cos(t5)*sin(s5), sin(t5)*sin(s5), cos(s5);
       cos(t6)*sin(s6), sin(t6)*sin(s6), cos(s6)];

    dm = acos(max(min(X*X',1),-1));
    [tb,td] = persistence_matrix(dm);
    bd_times(i,:) = [tb,td];
end

% Add a couple of points to make the graph prettier
bd_times = [pi/2,pi;
            bd_times(1:end-1,:);
            acos(-1/3), acos(-1/3)];

% I start from 40 to obtain only the left boundary. Other configurations
% have points inside of the region
% figure();
% hold on;

% plot(bd_times(:,1), bd_times(:,2), color);
% plot_reference_lines();

plot(2*sin(bd_times(:,1)/2), 2*sin(bd_times(:,2)/2), color);

tt=0:0.01:2;
plot(tt,tt,'red');
plot(sqrt(3/4)*tt,tt,'red');

tt=sqrt(2):0.01:2;
plot(tt,2*ones(size(tt)),'red');

% tt = sqrt(2):0.01:sqrt(8/3);
% tt = [tt, sqrt(8/3)];
% plot(tt, sqrt(4*tt.^2.*(3-tt.^2)./(4-tt.^2)), 'g.');

% hold off;
end