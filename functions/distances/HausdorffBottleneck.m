% Calculates the Hausdorff Distance between P and Q induced by the
% bottleneck distance.
% ---- Each row of P and Q is a single persistence diagram ----
function hd = HausdorffBottleneck(P,Q)
sP = size(P,1); sQ = size(Q,1);

% Precompute the persistence of all points
persP = P(:,2)-P(:,1);
persQ = Q(:,2)-Q(:,1);

dPQ = zeros(sP,1);
parfor p=1:sP
    % calculate the minimum bottleneck distance from points in P to Q
    d_inf = max( abs(P(p,:)-Q), [], 2);
    dPQ(p) = min(min(d_inf, 0.5*max(persP(p), persQ)));
end

dQP = zeros(sQ,1);
parfor q=1:sQ
    % calculate the minimum bottleneck distance from points in Q to P
    d_inf = max( abs(P-Q(q,:)), [], 2);
    dQP(q) = min(min(d_inf, 0.5*max(persP, persQ(q))));
end

% The Hausdorff distance is max(max_a d(a,A), max_b d(A,b))
hd = max([dPQ; dQP]);
