clear; close all; clc;
addpath('../../functions/persistence/')

% Parameters
b = 1;
% a = sqrt(2);
a = 2.5;

k=1;            % Dimension of homology to calculate
n=2*k+2;        % Number of points to use
nReps = 1e5;    % Number of snapshots to take

% Saving results
save_to_file = false;
results_file = sprintf('D41(Ellipse) a=%1.1f, b=%1.1f', a, b);
image_file = results_file;

% Set random seed
rng(1)

% Arrays to save results
confs = zeros(n,2,nReps);
dms = zeros(n,n,nReps);
bd_times = zeros(nReps,2);

for i=1:nReps
    % Sample uniformly from the interval [0,2*pi]
    T = 2*pi*rand(n,1);

    % Use the polar form of the ellipse
    X = [a*cos(T), b*sin(T)];

    % Calculate the distance matrix
    dm = L2_distance(X',X');

    % Store the sample
    confs(:,:,i) = X;
    dms(:,:,i) = dm;

    % Calculate homology
    [tb,td] = persistence_matrix(dm);
    bd_times(i,:) = [tb,td];
end

%% Persistence of configurations in a single half
t1 = 0:0.01:pi/2;
n1_r = length(t1);

tb_r = zeros(n1_r,n1_r);
td_r = zeros(n1_r,n1_r);

% figure();
hold on;

for i=1:n1_r
    t0 = t1(i);
    t2 = t0:0.01:pi/2;

    x1 = a*cos(t0)*ones(size(t2));
    x2 = a*cos(t0)*ones(size(t2));
    x3 = a*cos(t2);
    x4 = a*cos(t2);

    y1 =  sqrt(1-x1.^2/a^2);
    y2 = -sqrt(1-x2.^2/a^2);
    y3 = -sqrt(1-x3.^2/a^2);
    y4 =  sqrt(1-x4.^2/a^2);
    
    % x1.^2/a^2+y1.^2
    
    d12 = sqrt((x1-x2).^2+(y1-y2).^2);
    d23 = sqrt((x2-x3).^2+(y2-y3).^2);
    d34 = sqrt((x3-x4).^2+(y3-y4).^2);
    d41 = sqrt((x4-x1).^2+(y4-y1).^2);

    d13 = sqrt((x1-x3).^2+(y1-y3).^2);
    d24 = sqrt((x2-x4).^2+(y2-y4).^2);
    
    tb_r(i,i:end) = max([d12;d23;d34;d41], [], 1);
    td_r(i,i:end) = min(d13,d24);
    
    I = tb_r(i,:) < td_r(i,:);
    plot(tb_r(i,I), td_r(i,I), 'g.');
end
hold off;

%% Persistence of parallelograms
T1 = 0:0.01:pi/2;
T2 = 0:0.01:pi;

n1_p = length(T1);
n2_p = length(T2);

tb_p = zeros(n1_p,n2_p);
td_p = zeros(n1_p,n2_p);

% figure();
hold on;
for i=1:n1_p
    t1 = T1(i);
    
    I = t1<T2 & T2<pi-t1;
    t2 = T2(I);

    x1 = a*cos(t1)*ones(size(t2));
    x2 = a*cos(t2);
    x3 = a*cos(pi+t1)*ones(size(t2));
    x4 = a*cos(pi+t2);

    y1 =  sqrt(1-x1.^2/a^2);
    y2 =  sqrt(1-x2.^2/a^2);
    y3 = -sqrt(1-x3.^2/a^2);
    y4 = -sqrt(1-x4.^2/a^2);
    
    % x1.^2/a^2+y1.^2
    
    d12 = sqrt((x1-x2).^2+(y1-y2).^2);
    d23 = sqrt((x2-x3).^2+(y2-y3).^2);
    d34 = sqrt((x3-x4).^2+(y3-y4).^2);
    d41 = sqrt((x4-x1).^2+(y4-y1).^2);

    d13 = sqrt((x1-x3).^2+(y1-y3).^2);
    d24 = sqrt((x2-x4).^2+(y2-y4).^2);
    
    tb_p(i,I) = max([d12;d23;d34;d41], [], 1);
    td_p(i,I) = min(d13,d24);
    
    I = tb_p(i,:) < td_p(i,:);
    plot(tb_p(i,I), td_p(i,I), 'c.');
end
hold off;

%% Persistence of configurations close to a kite
% x0 = a/sqrt(a^2+1);
% t0 = acos(x0/a);

tt = 0:0.01:pi/6;
x1 = a*ones(size(tt));
x2 = a*cos(pi/2  +1*tt);
x3 =-a*ones(size(tt));
% x3 = a*cos(pi    +2*tt);
x4 = a*cos(3*pi/2-1*tt);

y1 = zeros(size(tt));
y2 = sin(pi/2  +1*tt);
% y3 = sin(pi    +2*tt);
y3 = zeros(size(tt));
y4 = sin(3*pi/2-1*tt);

d12 = sqrt((x1-x2).^2+(y1-y2).^2);
d23 = sqrt((x2-x3).^2+(y2-y3).^2);

d13 = sqrt((x1-x3).^2+(y1-y3).^2);
d24 = sqrt((x2-x4).^2+(y2-y4).^2);

tb = max(d12,d23);
td = min(d13,d24);

% figure();
hold on;

tt = 0:0.01:max(max(tb,td));

% plot(tt,ones(size(tt))*2,'red')
% plot(ones(size(tt))*sqrt(244)/10,tt,'red')

plot(tt,tt,'red');
plot(tb,td,'red');

hold off;

%% Plot the results
% figure();
hold on;

tt = 0:0.01:max(bd_times(:));
plot(tt,tt,'red');
plot(tt,sqrt(2)*tt,'red');

plot(bd_times(:,1), bd_times(:,2), 'b.');

title(sprintf('D41(Ellipse) a=%1.1f, b=%1.1f', a, b));
hold off;

%% Persistence of rectangles
x0 = a/sqrt(a^2+1);
y0 = x0;

tt = 0:0.01:(a-x0);
x1 =  x0+tt;
x2 = -x0-tt;
x3 = -x0-tt;
x4 =  x0+tt;

y1 =  sqrt(1-x1.^2/a^2);
y2 =  sqrt(1-x2.^2/a^2);
y3 = -sqrt(1-x3.^2/a^2);
y4 = -sqrt(1-x4.^2/a^2);

d12 = sqrt((x1-x2).^2+(y1-y2).^2);
d23 = sqrt((x2-x3).^2+(y2-y3).^2);

d13 = sqrt((x1-x3).^2+(y1-y3).^2);
d24 = sqrt((x2-x4).^2+(y2-y4).^2);

tb = max(d12,d23);
td = min(d13,d24);

% figure();
hold on;

tt = 0:0.01:max(max(tb,td));
plot(tt,tt,'red');
plot(tb,td,'red');

hold off;

%% Saving results
% Save the results
if save_to_file
    saveas(gcf, sprintf('./results/%s.png', image_file)); %#ok<UNRCH>
    save(sprintf('./results/%s.mat',results_file));
end
