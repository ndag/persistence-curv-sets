% Compute 0-dimensional VR-persistent homology
% INPUT: N-by-N distance matrix dm
% OUTPUT: (N-1)-dim'l row vector with the death times of the finite
% intervals.
% I don't output birth times because they are all 0, and I don't output the
% infinite interval because every finite set has one.
function dgm = dgm_0_um(dm)
N = size(dm,1);
dmc = dm;
clusters = num2cell(1:N);
dgm = zeros(1,N-1);
nInts = 0;

% We'll keep going as long as there is more than 1 cluster
while numel(clusters)>1
    n = size(dmc,1);
    
    % Find the smallest distance
    r = min(squareform(dmc));
    
    % Find clusters
    I_sep = 1:n;
    cl_new = cell(0,0);
    idx_cl = 1;
    while ~isempty(I_sep)
        % Find all points that are at distance r from point i
        i = min(I_sep);
        cluster = [i, find(dmc(i,:)==r)];
        
        % Add the cluster to the list
        cl_new{idx_cl} = cluster;
        idx_cl = idx_cl+1;
        
        % Remove the clustered elements from I_sep
        I_sep = setdiff(I_sep, cluster);
    end
    
    % Update the global cluster list
    cl_old = clusters;
    nClNew = numel(cl_new);
    nClOld = numel(cl_old);
    
    clusters = cell(1,nClNew);
    for idx_cl=1:nClNew
        cl = cl_new{idx_cl};
        clusters{idx_cl} = [cl_old{cl}];
    end
    
    % Compute a distance matrix between the new clusters
    dmc_old = dmc;
    dmc = zeros(nClNew, nClNew);
    for i=1:nClNew
        cl1 = cl_new{i};
        for j=i+1:nClNew
            cl2 = cl_new{j};
            
            % The distance between clusters 1 and 2 is the minimum distance
            % d(x1, x2) where xi is in cluster i
            dm_cls = dmc_old(cl1, cl2);
            dmc(i,j) = min(dm_cls(:));
            dmc(j,i) = dmc(i,j);
        end
    end
    
    % Add intervals with death time r to dgm_0
    nIntsNew = nClOld-nClNew;
    dgm(nInts+(1:nIntsNew)) = r*ones(1,nIntsNew);
    nInts = nInts+nIntsNew;
end

end