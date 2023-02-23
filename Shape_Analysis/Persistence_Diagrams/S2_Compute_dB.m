% Compute Bottleneck distances
% To compute dBs in parallel, I need to write a single parfor loop because
% nested loops are troublesome. For that, I'll compute squareform(dBs) and
% reassemble the square matrix dBs when I need it.
% First, I need to write the list of indices.
I = 1:nShapes;
[Ix, Iy] = meshgrid(I,I);
I = [Ix(:), Iy(:)];
I = I( I(:,1) < I(:,2), :);
Ix = I(:,1);
Iy = I(:,2);

% In I, the y coordinate advances before the x coordinate, i.e.
% I = [1,2; 1,3; 1,4; ...; 2,3; 2,4; ...]
% To verify that squareform assembles the matrix as we want, change nShapes
% to 5 (to make the example small), recalculate I as above and run the next
% two lines:
% row = 1:nchoosek(nShapes,2); M = squareform(row);
% for i=1:size(I,1); disp(M(I(i,1), I(i,2))); end;
% This should display the numbers 1, ..., 10 in order.

% Now we compute squareform(dBs)
dBs = zeros(maxHomDim+1, nchoosek(nShapes,2));
parfor idx=1:numel(Ix)
    t = getCurrentTask();
    
    ix=Ix(idx);
    iy=Iy(idx);
    
    for d=1:N2
        PD1 = PD_all{ix,d}; %#ok<*PFBNS>
        PD2 = PD_all{iy,d};
        
        % Truncate persistence diagrams
        thresh = min(T_all(ix), T_all(iy));
        PD1 = TruncatePD(PD1, thresh, false);
        PD2 = TruncatePD(PD2, thresh, false);
        
        tic
        dB = BottleneckDistancePar(PD1, PD2, t.ID);
        dBs(d,idx) = dB;
        dt = toc;
        fprintf('(%i,%i)\n dB = %0.5f\n Time: %0.5f\n\n ', idx,size(I,1), dB, dt)
    end
end

% The final distance matrix C is the maximum of dB over all dimensions
C = max(max(dBs),0);

% Save results
for d=1:maxHomDim+1
    mFile2 = matfile(sprintf('%s/Distances/dB_%i.mat', results_folder, d-1), 'Writable', true);
    mFile2.dB = dBs(d,:);
end
mFile2 = matfile(sprintf('%s/Distances/dB_max.mat', results_folder), 'Writable', true);
mFile2.dB = C;
