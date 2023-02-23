function G = buildGraph_Sumner(surface)
nT = length(surface.f.v);
nV = length(surface.v);

G = sparse(nV,nV);

for k=1:nT
    ak = surface.f.v(k,1);
    bk = surface.f.v(k,2);
    ck = surface.f.v(k,3);
    
    G(ak,bk) = 1;
    G(ak,ck) = 1;
    G(bk,ck) = 1;
end

% Add weights edge by edge
[ii,jj,~] = find(G);
nE = length(ii);
for k=1:nE
    P1 = surface.v(ii(k), :);
    P2 = surface.v(jj(k), :);
    G(ii(k), jj(k)) = L2_distance(P1', P2');
end

% symmetrize
G = max(G,G');
