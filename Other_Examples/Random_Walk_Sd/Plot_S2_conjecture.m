clear; close all; clc;
addpath('../../functions/persistence/');
addpath('../../functions/functions_S2/');

% Load the known persistence to compare it
file_name = 'results/D62_S2_Random_Walk.mat';
myMatfile = matfile(file_name);
bd_times_n = myMatfile.bd_times_n;
bd_times_0 = myMatfile.bd_times_0;

save_to_file = false;
results_file = 'D62_S2_Conjecture';
image_file = 'D62_S2_Conjecture';

% Parameters for the simulation
n = 6;              % Number of points in each configuration
d = 2;              % Dimension of the sphere

% Number of configurations to sample
numH1 = 150;
numH2 = 100;
numT = 100;

% Lists to store the results
confs = zeros(n,d+1,numT,numH1,numH2);
dms = zeros(n,n,numT,numH1,numH2);
bd_times = zeros(numT,2,numH1,numH2);

% Create triangle configurations
minH = 0;
HH = minH:((pi/2-minH)/(numH1-1)):pi/2;
if numT>1
    TT = 0:((pi/3)/(numT-1)):pi/3;
else
    % TT = pi/3;
    TT = pi/3;
end

%% Perform the experiment
% HH = 0.2;
% TT = pi/6;

for i=1:numT
    disp(i)
    for j=1:numH1
        T = TT(i);
        H = HH(j);
        
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
        HH2 = -H:((pi/2+H)/(numH2-1)):pi/2;
        % HH2 = 0:((pi/2)/(numH2-1)):pi/2;
        % HH2 = -H:H/(numH2-1):0;
        % HH2 = -minH;
        if ~isempty(HH2)
        for k=1:numH2
            H2 = HH2(k);
            s4 = pi/2+H2;
            s5 = pi/2+H2;
            s6 = pi/2+H2;

            % Construct the configuration
            X = [cos(t1)*sin(s1), sin(t1)*sin(s1), cos(s1);
                 cos(t2)*sin(s2), sin(t2)*sin(s2), cos(s2);
                 cos(t3)*sin(s3), sin(t3)*sin(s3), cos(s3);
                 cos(t4)*sin(s4), sin(t4)*sin(s4), cos(s4);
                 cos(t5)*sin(s5), sin(t5)*sin(s5), cos(s5);
                 cos(t6)*sin(s6), sin(t6)*sin(s6), cos(s6)];
        
            % Calculate the distance matrix
            dm = L2_distance(X',X');
            dm = 2*asin(dm/2);

            % Calculate persistence
            [tb,td] = persistence_matrix(dm);

            % Store the results
            confs(:,:,i,j,k) = X;
            dms(:,:,i,j,k) = dm;
            bd_times(i,:,j,k) = [tb,td];
        end
        end % if ~isempty(HH2)
    end
end

B0 = zeros(numT*numH1*numH2, 2);
B0(:,1) = reshape(bd_times(:,1,:,:), [], 1);
B0(:,2) = reshape(bd_times(:,2,:,:), [], 1);

%% Plot the persistence
figure();
hold on;

% Plot the persistence we just found
for z=1:size(bd_times,3)
    for w=1:size(bd_times,4)
        plot(2*sin(bd_times(:,1,z,w)/2), 2*sin(bd_times(:,2,z,w)/2), 'b.');
    end
end

% Plot the known persistence to compare
plot(2*sin(bd_times_n(:,1)/2), 2*sin(bd_times_n(:,2)/2), 'm.');
plot(2*sin(bd_times_0(:,1)/2), 2*sin(bd_times_0(:,2)/2), 'g.');
plot_boundary_configurations();

% Euclidean
ylim([0,2.25])
text(0.1,2.20, 'Green = Uniform sample; Magenta = Random Walk');
text(0.1,2.10, 'Blue = Obtained from two parallel equilateral triangles.');
title('$$D_{6,2}^\mathrm{VR}(S^2_E)$$','interpreter','latex')

% Format
set(gca,'FontSize',16)

hold off;

%% Save results
if save_to_file
    saveas(gcf, sprintf('results/%s.png', image_file));
end
