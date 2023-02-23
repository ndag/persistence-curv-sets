% Compute persistence sets and histograms
disp('Started PS0')
parfor idx=1:nShapes
    tic
    
    name = names{idx};
    
    file_name = sprintf('%s/%s.mat', input_folder, name);
    dm = load(file_name);
    
    if size(dm.dm,1)>1
        N = size(dm.dm,1);
        dm = squareform(dm.dm);
    else
        dm = dm.dm;
        N = size(squareform(dm),1);
    end
    
    % The principal persistence set D_{2,0} equals the distribution of
    % distances of a metric space
    % Recall that we can also choose x1=x2, so D_{2,0} also contains N
    % copies of the 0 diagram
    dm = [zeros(1,N), dm];
    bd_times = [zeros(length(dm),1), sort(dm)'];
    
    dt = toc;
    fprintf('(%i/%i) %s (calc): %0.2f\n ', idx,nShapes,name,dt)
    
    % Save results
    tic
    file = matfile(sprintf('%s/D_%i_%i/%s.mat', results_folder, n,k, name), 'Writable', true);
    file.bd_times = bd_times;
    
    dt = toc;
    fprintf('(%i/%i) %s (save): %0.2f\n ', idx,nShapes,name,dt)
end
