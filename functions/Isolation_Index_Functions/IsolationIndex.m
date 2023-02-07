function [alpha, A, B] = IsolationIndex(dm, I, divide)
    if nargin<3
        divide = true;
    end
    
    n = size(dm,1);
    J = setdiff(1:n, I);
    
    alpha = Inf;
    A = 0;
    B = 0;
    for i1=1:length(I)
        for i2=i1:length(I)
            for j1=1:length(J)
                for j2=j1:length(J)
                    a1 = I(i1);
                    a2 = I(i2);
                    b1 = J(j1);
                    b2 = J(j2);
                    
                    d1 = dm(a1,b1)+dm(a2,b2);
                    d2 = dm(a1,b2)+dm(a2,b1);
                    d3 = dm(a1,a2)+dm(b1,b2);
                    
                    Beta = max([d1,d2,d3])-d3;
                    
                    if Beta < alpha
                        alpha = Beta;
                        A = [a1,a2];
                        B = [b1,b2];
                    end
                end
            end
        end
    end
    
    % Sometimes I want to make the division by 2 until the end of a long
    % calculation. If divide==true, I want the division now.
    if divide
        alpha = 0.5*alpha;
    end
end