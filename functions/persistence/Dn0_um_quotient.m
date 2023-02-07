% Computes the collection {D_{n,0}(U/~R)}_n of the quotient ultrametric space
% U/~R, where a ~R b iff d(a,b)<=R. In terms of the D_{n,0}, this amounts
% to eliminating the intervals with death time smaller than R
function Dn0_q = Dn0_um_quotient(Dn0,R)
% Compute the total number of diagrams
nDiags = sum(cellfun(@(D) size(D,1), Dn0));
Dn0_new = cell(nDiags,1);

idx=1;
% For each number of points
for n=1:numel(Dn0)
    Dn = Dn0{n};
    
    % For each diagram in Dn, only keep the entries above R
    L = size(Dn,1);
    for i=1:L
        D = Dn(i,:);
        Dn0_new{idx} = D(D>R);
        idx = idx+1;
    end
end

% Sort by size
sizes = cellfun(@(D) numel(D), Dn0_new)+1;
Dn0_q = cell(max(sizes), 1);
for idx=1:nDiags
    n = sizes(idx);
    if n>0
        Dn0_q{n} = [Dn0_q{n}; Dn0_new{idx}];
    end
end

% Remove duplicates
Dn0_q = cellfun(@(M) unique(M,'rows'), Dn0_q, 'UniformOutput', false);
end