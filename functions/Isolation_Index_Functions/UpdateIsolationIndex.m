% alpha_0 is the isolation index of the partial dm-split A,B. This
% function calculates the isolation index of the partial split
% A \cup {x_A}, B by finding the minimum of Beta_{a,x_A},{b1,b2} for all
% a \in A and b1, b2 \in B. Then the new alpha is the minimum of that
% quantity and alpha_0.
% Notice that we don't calculate Beta_{a1,a2},{b1,b2} with a_2 \neq x
% because that was already calculated to obtain alpha_0.
% IMPORTANT: This function actually calculates twice the isolation index. I
% don't divide by 2 to avoid unnecessarily dividing or multiplying by 2
% during each iteration.
function alpha = UpdateIsolationIndex(x_A, A, B, alpha_0, dm)
    alpha = alpha_0;
    a2 = x_A;
    
%     for i1=1:length(A)+1
    for i1=1:length(A)
        for j1=1:length(B)
            for j2=j1:length(B)
%                 if i1<length(A)
%                     a1 = A(i1);
%                 else
%                     a1 = a2;
%                 end
                a1 = A(i1);
                b1 = B(j1);
                b2 = B(j2);

                d1 = dm(a1,b1)+dm(a2,b2);
                d2 = dm(a1,b2)+dm(a2,b1);
                d3 = dm(a1,a2)+dm(b1,b2);

                Beta = max([d1,d2,d3])-d3;

                alpha = min(alpha, Beta);
            end
        end
    end
end