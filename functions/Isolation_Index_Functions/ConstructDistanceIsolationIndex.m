function dm_a = ConstructDistanceIsolationIndex(alpha, A, T, gamma)
    if nargin==3
        gamma=0;
    end
    
    n = length(alpha);
    M_alpha = alpha+alpha'-2*diag(alpha);
    MA = zeros(size(A));
    
    if T==1
        for i=1:n
            i1 = n-mod(n-(i+1),n);      % i+1 mod n with representative in 1:n
            MA = MA+A(i,i1)*SplitMetric([i,i1],n);
        end
    elseif T==2
        % We assume 1:5=[x,y,u,v,w]
        MA = A(1,3)*SplitMetric([1,3],n)...
            + A(1,4)*SplitMetric([1,4],n)...
            + A(2,3)*SplitMetric([2,3],n)...
            + A(2,4)*SplitMetric([2,4],n);
        
        D2 = ones(n,n) - eye(n,n);
        D2(1,2) = 2;
        D2(2,1) = 2;
        D2(3,4) = 2;
        D2(4,3) = 2;
        D2(3,5) = 2;
        D2(5,3) = 2;
        D2(4,5) = 2;
        D2(5,4) = 2;
        
        MA = MA+gamma*D2;
    elseif T==3
        % We assume 1:5=[x,y,u,v,w]
        MA = A(1,3)*SplitMetric([1,3],n)...
            + A(1,4)*SplitMetric([1,4],n)...
            + A(2,4)*SplitMetric([2,4],n)...
            + A(2,5)*SplitMetric([2,5],n);
        
        D2 = ones(n,n) - eye(n,n);
        D2(1,2) = 2;
        D2(2,1) = 2;
        D2(3,4) = 2;
        D2(4,3) = 2;
        D2(3,5) = 2;
        D2(5,3) = 2;
        D2(4,5) = 2;
        D2(5,4) = 2;
        
        MA = MA+gamma*D2;
    end
    
    dm_a = M_alpha+MA;
end