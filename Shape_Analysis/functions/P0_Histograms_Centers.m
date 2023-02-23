
%% Histogram for k=0
nBins = nBins0;

Cs = linspace(0,1,nBins);
Cs = [zeros(nBins,1), Cs'];

% Compute cost matrix
C = d_bottleneck(Cs, Cs);

% Save parameters to disk
save(sprintf('%s/Hist_Dn0.mat', param_folder),...
    'nBins','Cs','C',...
    '-v7.3');

%% Histograms for k>1
% Number of squares to distribute the points into
nBins = nBins1;
nBinsDim = ceil(sqrt(2*nBins));
disp(nBinsDim^2)

% Maximum values in each dimension
Mx = 1;
My = 1;

% Centers
Cx = linspace(0,Mx, nBinsDim);
Cy = linspace(0,My, nBinsDim);

% 2d vector of centers
[xx,yy] = meshgrid(Cx,Cy);

% Transpose to have the same order as HistBD
xx = xx'; yy = yy';

% Expand into vector of centers and compute the cost matrix (ie distance matrix)
Cs = [xx(:), yy(:)];

% Only keep points above the diagonal and 0
Ic = (Cs(:,1) < Cs(:,2)) | sum(Cs,2)==0;
disp(sum(Ic))
Cs = Cs(Ic,:);

% Compute cost matrix
C = d_bottleneck(Cs, Cs);

% In Cs, x increases, then y does
% Matlab expands matrix by transversing rows then columns. Each column
% HistBD(:,j,idx,idy) has y fixed while x increases, so we need to define
% Cs as we did so that the variables are synced when expanding M(:)

% Save parameters to disk
save(sprintf('%s/Hist_%i.mat', param_folder, size(C,1)),...
    'nBins', 'nBinsDim', 'Mx', 'My',...
    'Cx','Cy','Cs','Ic','C',...
    '-v7.3');
