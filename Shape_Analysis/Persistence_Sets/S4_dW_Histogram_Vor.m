% Load cost matrix
C = Params.C;

% Load histograms
Hists = zeros(nShapes,size(C,1));
tic
for idx=1:nShapes
    % Format file name
    name = names{idx};
    
    % Load histogram. I don't need to filter out any points because I
    % constructed the centers manually (instead of using hist3)
    mFile = matfile(sprintf('%s/D_%i_%i/%s.mat', results_folder,n,k, name));
    Hists(idx,:) = mFile.Hist_Vor;
end
dt=toc;
fprintf('Loading: %0.5f\n ', dt)

% % To resume after calculating on a subset
% dw_file = matfile(sprintf('%s/Distances/dW_%i_%i.mat', results_folder,n,k), 'Writable', true);
% dWs = dw_file.dW;

% To compute dWs in parallel, I need to write a single parfor loop because
% nested loops are troublesome. For that, I'll compute squareform(dWs) and
% I'll reassemble the square matrix dWs after I'm done.
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
dWs = zeros(1,nchoosek(nShapes,2));
parfor idx=1:numel(Ix)
    ix=Ix(idx);
    iy=Iy(idx);
    
    hist_1 = Hists(ix,:);
    hist_2 = Hists(iy,:);
    
    tic;
    try
        [dW, ~] = emd_mex(hist_1,hist_2, C);
    catch EM
        fprintf('%s (idx,ix,iy) = (%i,%i,%i)\n Will set dW as infinity', EM.message, idx,ix,iy);
        dW = inf;
    end
    dt = toc;
    dWs(idx) = dW;
    fprintf('(%i,%i)\n dW = %0.5f\n Time: %0.5f\n\n ', idx,size(I,1), dW, dt)
end

% Save results
if useShortlist
    dW_file = matfile(sprintf('%s/Distances/dW_%i_%i_(%i).mat', results_folder,n,k,nShapes), 'Writable', true);
else
    dW_file = matfile(sprintf('%s/Distances/dW_%i_%i.mat', results_folder,n,k), 'Writable', true);
end

dW_file.dW = dWs;

clear dW_file
