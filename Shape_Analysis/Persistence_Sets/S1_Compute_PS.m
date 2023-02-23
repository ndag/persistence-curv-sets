% Compute persistence sets and histograms
parfor idx=1:nShapes
    tic
    
    % Lists to save results
    confs = zeros(nReps,n);
    bd_times = zeros(nReps,2);
    
    name = names{idx};

    file_name = sprintf('%s/%s.mat', input_folder, name);
    dm = load(file_name);
    dm0 = squareform(dm.dm);
    N = size(dm0,1);

    % Calculate the persistence set
    for i=1:nReps
        I = randi(N,n,1);
        dm = dm0(I,I);
        [tb,td] = persistence_matrix(dm);
        bd_times(i,:) = [tb,td];
        confs(i,:) = I';
    end
    
    dt = toc;
    fprintf('(%i/%i) %s (calc): %0.2f\n ', idx,nShapes,name,dt)
    
    % Save results
    tic
    file = matfile(sprintf('%s/D_%i_%i/%s.mat', results_folder, n,k, name), 'Writable', true);
    file.confs = confs;
    file.bd_times = bd_times;
    
    dt = toc;
    fprintf('(%i/%i) %s (save): %0.2f\n ', idx,nShapes,name,dt)
end
