% Take a sample of N points from each space
for idx=1:nShapes
    tic
    name = names{idx};
    
    file_name = sprintf('%s/%s.mat', input_folder, name);
    mFile = matfile(file_name, 'Writable', true);
    dm = squareform(mFile.dm);
    
    % Subsample and renormalize
    I = px_fps(dm,'metric', nL,'n');
    dmX = dm(I,I);
    dmX = dmX/max(dmX(:));
    
    % Compute the threshold for the persistence diagrams
    thresh = SelectThreshold(dmX, s_edges);
    
    % Save the results
    mFile.I = I;
    mFile.dmX = squareform(dmX);
    mFile.thresh = thresh;
    
    dt = toc;
    fprintf('(%i/%i) done: %0.5f\n ', idx,nShapes,dt);
end
