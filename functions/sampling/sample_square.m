function [X, dm] = sample_square(N, p, R1, R2)
    % We construct a rectangle [0,2*R1]*[0,2*R2].
    % Right now, only the L_2 and L_inf metrics are programmed as options.
    if nargin<2
        p=inf;
    end
    if nargin<4
        R1=pi;
        R2=pi;
    end

    X = rand(N,2);
    X(:,1) = 2*R1.*X(:,1);
    X(:,2) = 2*R2.*X(:,2);
    
    if ~isinf(p)
        dm = L2_distance(X',X');
    else
        d1 = abs(X(:,1)-(X(:,1))');
        d2 = abs(X(:,2)-(X(:,2))');
        dm = max(d1,d2);
    end
end


%     X = zeros(N,3);
%     
%     theta_1 = 2*pi.*rand(N,1);
%     theta_2 = 2*pi.*rand(N,1);
%     
%     X(:,1) = (R + r*cos(theta_2)).*cos(theta_1);
%     X(:,2) = (R + r*cos(theta_2)).*sin(theta_1);
%     X(:,3) = r*sin(theta_2);
%     
%     % We'll just do Euclidean distance for now
%     dm = L2_distance(X',X');