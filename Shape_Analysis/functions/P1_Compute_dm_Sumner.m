folders = dir('models');
folders = folders([folders.isdir]);
folders = {folders(3:end).name};
% Remove the folder with processed matrices
folders = folders( ~strcmp(folders, 'processed') );
nFolders = numel(folders);

for idx=1:nFolders
    % Read the names of .obj files in the folder
    folder = folders{idx};
    files_dir = dir(sprintf('models/%s/*.obj', folder));
    files_dir = {files_dir.name};
    
    % Remove extension
    files_dir = cellfun(@(str) str(1:end-4), files_dir,...
        'UniformOutput', false);
    
    % Calculate distance matrices
    nFiles = numel(files_dir);
    
    fprintf('Start: %s\n', folder);
    for idf=1:nFiles
        tic
        name = files_dir{idf};
        shape = readObj(sprintf('models/%s/%s.obj', folder, name));
        G = buildGraph_Sumner(shape);
        
        % find largest connected component of the graph (some points can be
        % disconnected)
        [aa, bb] = conncomp(graph(G),'OutputForm','cell');

        fprintf('There are %i points in the largest connected component\n', bb(1))
        fprintf('and %i in the full shape\n', size(G,1));
        Jk = aa{1};
        
        % get the point cloud
        P = shape.v;
        
        % retain only points in largest component of the graph
        P = P(Jk,:);
        G = G(Jk,Jk);

        % apply FPS
        I = px_fps(P','vector',NFPS,'n');
        
        % compute the geodesic distance matrix
        dm = distances(graph(G),I,I);
        dm = squareform(dm);
        
        % Normalize distances
        dm = dm/max(dm(:));
        
        dt = toc;
        fprintf('(%i,%i)/(%i,%i) %s (calc): %0.2f\n ', idx,idf,nFolders,nFiles,name,dt)
        
        % Save results
        tic
        mFile = matfile(sprintf('%s/%s.mat', input_folder,name), 'Writable', true);
        mFile.dm = dm;
        dt = toc;
        fprintf('(%i,%i)/(%i,%i) %s (save): %0.2f\n\n', idx,idf,nFolders,nFiles,name,dt)
    end
end
