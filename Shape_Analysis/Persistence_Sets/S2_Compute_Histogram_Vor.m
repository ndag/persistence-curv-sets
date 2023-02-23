% Compute Voronoi histograms
% Load the histogram centers
Cs = Params.Cs;

% I add each point to the bin whose center is the closest
parfor idx=1:nShapes
    tic
    
    % Placeholder for Histograms
    Hist_Vor = zeros(size(Cs,1),1);
    
    % Load the persistence set we already calculated
    name = names{idx};
    
    mFile = matfile(sprintf('%s/D_%i_%i/%s.mat', results_folder,n,k, name), 'Writable', true);
    bd_times = mFile.bd_times;
    
    % Compute the distances from each point in bd_times to the centers
    nReps = size(bd_times,1);
    for i=1:nReps
        dgm = bd_times(i,:);
        dB = d_bottleneck(dgm, Cs);
        
        % If there are several minima, this code chooses the one with the
        % lowest x and y coordinates.
        [~,I] = min(dB);
        Hist_Vor(I) = Hist_Vor(I)+1;
        
        % If there are several minima, this code gives more control over
        % which one to use.
        % [~,y] = find(dB==min(dB))
        % Hist_Vor(y(1)) = Hist_Vor(y(1))+1;
    end
    
    dt = toc;
    fprintf('(%i/%i) %s (calc): %0.2f\n ', idx,nShapes,name,dt)
    
    % Save the histogram back to the same file
    tic
    mFile.Hist_Vor = Hist_Vor;
    dt = toc;
    fprintf('(%i/%i) %s (save): %0.2f\n ', idx,nShapes,name,dt)
end
