% Load histograms
bd_cell = cell(nShapes,1);

tic
for idx=1:nShapes
    % Format file name
    name = names{idx};
    
    % Load histogram
    mFile = matfile(sprintf('%s/D_%i_%i/%s.mat', results_folder,n,k, name));
    
    % Removing duplicate diagrams greatly speeds up the computation of the
    % Hausdorff distance
    bd_cell{idx} = unique(mFile.bd_times, 'rows');
end
dt=toc;
fprintf('Loading: %0.5f\n ', dt)

% % To resume after calculating a subset
% dH_file = matfile(sprintf('%s/dW/dH_%i_%i.mat', results_folder,n,k), 'Writable', true);
% dHs = dH_file.dH;

% To compute dHs in parallel, I need to write a single parfor loop because
% nested loops are troublesome. For that, I'll compute squareform(dHs) and
% I'll reassemble the square matrix dHs after I'm done.
% First, I need to write the list of indices.
I = 1:nShapes;
[Ix, Iy] = meshgrid(I,I);
I = [Ix(:), Iy(:)];
I = I( I(:,1) < I(:,2), :);
Ix = I(:,1);
Iy = I(:,2);
N = numel(Ix);

% In I, the y coordinate advances before the x coordinate, i.e.
% I = [1,2; 1,3; 1,4; ...; 2,3; 2,4; ...]
% To verify that squareform assembles the matrix as we want, change nShapes
% to 5 (to make the example small), recalculate I as above and run the next
% two lines:
% row = 1:nchoosek(nShapes,2); M = squareform(row);
% for i=1:size(I,1); disp(M(I(i,1), I(i,2))); end;
% This should display the numbers 1, ..., 10 in order.

% Now we compute squareform(dHs)
dHs = zeros(1,nchoosek(nShapes,2));
parfor idx=1:numel(Ix)
    ix=Ix(idx);
    iy=Iy(idx);
    
    bd_1 = bd_cell{ix,1};
    bd_2 = bd_cell{iy,1};
    
    tic;
    dH = HausdorffBottleneck(bd_1,bd_2);
    dt = toc;
    dHs(idx) = dH;
    fprintf('(%i,%i)\n dH = %0.5f\n Time: %0.5f\n\n ', idx,size(I,1), dH, dt)
end

% Removed used variables
clear bd_cell

% Save results
if useShortlist
    dH_file = matfile(sprintf('%s/Distances/dH_%i_%i_(%i).mat', results_folder,n,k,nShapes), 'Writable', true);
else
    dH_file = matfile(sprintf('%s/Distances/dH_%i_%i.mat', results_folder,n,k), 'Writable', true);
end

dH_file.dH = dHs;

clear dH_file
