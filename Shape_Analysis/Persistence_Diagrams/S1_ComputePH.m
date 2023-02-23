N2 = maxHomDim+1;

% Lists to save persistence diagrams
PD_all = cell(nShapes,maxHomDim+1);
T_all = zeros(nShapes,1);

% Load distance matrices
parfor idx=1:nShapes
    t = getCurrentTask();
    % t.ID=1;
    tic
    
    name = names{idx};

    % Parameters for the metric space
    file_name = sprintf('%s/%s.mat', input_folder, name);
    mFile = matfile(file_name, 'Writable', true);
    dmX = squareform(mFile.dmX);
    T_all(idx) = mFile.thresh;
    
    % Compute PDs
    PDs = RipsFiltrationDMPar(dmX, maxHomDim, t.ID, '', T_all(idx));
    
    % Truncate and simplify diagrams
    for d=1:N2
        PD = PDs{d};
        PD = TruncatePD(PD, mFile.thresh, false);
        
        PD_all{idx,d} = PD;
    end
end

% Save PDs to disk
for idx=1:nShapes
    name = names{idx};
    
    file_name = sprintf('%s/%s.mat', input_folder, name);
    mFile = matfile(file_name, 'Writable', true);
    
    mFile.PDs = PD_all(idx,:);
end
