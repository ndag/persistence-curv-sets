% Calculates the Hausdorff Distance between P and Q induced by the
% bottleneck distance, assuming that both P and Q line on the line
% tb=constant
% This code accepts two input formats:
% 1. P and Q have shape [N,2], [M,2]. Since P(i,1)==Q(j,1) for any i and j,
%    we will retain only P(:,2) and Q(:,2).
% 2. P and Q have shape [N,1], [M,1]. The tb coordinate does not affect
% these calculations.
function hd = HausdorffBottleneckLine(P,Q)
if size(P,2) > 1
    P = P(:,2);
    Q = Q(:,2);
end

% Remove duplicate points. The unique function also sorts the arrays
P = unique(P);
Q = unique(Q);
sP = length(P); sQ = length(Q);

% Compute the minimum bottleneck distance for each point in P
dPQ = zeros(sP,1);

ip = 1;
iq = 1;
% I don't reset iq to 1 in each iteration because both P and Q are sorted.
while ip<=sP
    p = P(ip);
    
    % Find the first point in Q larger than p
    while iq<=sQ && Q(iq) <= p
        iq=iq+1;
    end
    
    if iq==1
        % All points in Q are larger than p
        d1 = inf;
        
        q2 = Q(iq);
        d2 = min(q2-p, q2/2);   % pers(p)<=pers(q2)=q2
    elseif iq==sQ+1
        % All points in Q are smaller than or equal to p
        q1 = Q(iq-1);
        d1 = min(p-q1, p/2);    % pers(q1)<=pers(p)=p
        
        d2 = inf;
    else
        % In non-edge cases
        q1 = Q(iq-1);
        q2 = Q(iq);

        % Calculate bottleneck distances
        d1 = min(p-q1, p/2);    % pers(q1)<=pers(p)=p
        d2 = min(q2-p, q2/2);   % pers(p)<=pers(q2)=q2
    end
    
    % Retain the minimum
    dPQ(ip) = min(d1,d2);
    ip=ip+1;
end

% Compute the minimum bottleneck distance for each point in Q
dQP = zeros(sQ,1);

ip = 1;
iq = 1;
% I don't reset ip to 1 in each iteration because both P and Q are sorted.
while iq<=sQ
    q = Q(iq);
    
    % Find the first point in P larger than q
    while ip<=sP && P(ip) <= q
        ip=ip+1;
    end
    
    if ip==1
        % All points in P are larger than q
        d1 = inf;
        
        p2 = P(ip);
        d2 = min(p2-q, p2/2);   % pers(q)<=pers(p2)=p2
    elseif ip==sP+1
        % All points in P are smaller than or equal to q
        p1 = P(ip-1);
        d1 = min(q-p1, q/2);    % pers(p1)<=pers(q)=q
        
        d2 = inf;
    else
        p1 = P(ip-1);
        p2 = P(ip);

        % Calculate bottleneck distances
        d1 = min(q-p1, q/2);    % pers(p1)<=pers(q)=q
        d2 = min(p2-q, p2/2);   % pers(q)<=pers(p2)=p2
    end
    
    % Retain the minimum
    dQP(iq) = min(d1,d2);
    iq=iq+1;
end

% The Hausdorff distance is max(max_p d(p,Q), max_q d(P,q))
hd = max([dPQ; dQP]);
