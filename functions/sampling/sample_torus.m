function [X, dm] = sample_torus(N, p, R1, R2)
    % We construct a torus modeled as the square [0,2*R1]x[0,2*R2] with the
    % boundaries identified in the usual way. This is homeomorphic to the
    % product (R1'*S^1)x(R2'*S^1), where Ri'=Ri/pi.
    % The metric is given by the function d_torus.
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
    
    dm = d_torus(X,X, p,R1,R2);
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