function [X,dm,T,RR] = ReorderToAcuteTriangle(X,dm, dmIsEuclidean)
    if nargin<3
        dmIsEuclidean = false;
    end
    n = size(X,1);
    
    % Get the persistence of the matrix
    [tb,~] = persistence_matrix(dm);
    
    % Reorder points so that d_12 = tb
    [row, col] = find(dm==tb);
    if row(1) < col(1)
        i1 = row(1); i2 = col(1);
    else
        i1 = col(1); i2 = row(1);
    end
    X([1,i1],:) = X([i1,1],:);
    dm([1,i1],:) = dm([i1,1],:);
    dm(:,[1,i1]) = dm(:,[i1,1]);
    
    X([2,i2],:) = X([i2,2],:);
    dm([2,i2],:) = dm([i2,2],:);
    dm(:,[2,i2]) = dm(:,[i2,2]);
    
    % Reorder points so that v_d(x_i) = x_i+3
    for i1=1:n/2
        i2 = find(dm(:,i1)==max(dm(:,i1)));

        X([i2,i1+3],:) = X([i1+3,i2],:);
        dm([i2,i1+3],:) = dm([i1+3,i2],:);
        dm(:,[i2,i1+3]) = dm(:,[i1+3,i2]);
    end

    % Find the triangle xi,xj,xk formed by non-opposite points with largest
    % circumradius
    TT = sort(setprod([1,4],[2,5],[3,6]),2);
    RR = zeros(size(TT,1),1);
    
    % Calculate the circumradii of the acute triangles
    if dmIsEuclidean
        dm_E = dm;
    else
        dm_E = 2*sin(dm/2);
    end
    
    for i=1:size(TT,1)
        T = TT(i,:);
        % Sides of the triangle
        a = dm_E(T(2),T(3)); b = dm_E(T(3),T(1)); c = dm_E(T(1),T(2));
        % Cosines given by the cosine law
        cos_T = [(-a^2+b^2+c^2)/(2*b*c),...
                 (+a^2-b^2+c^2)/(2*a*c),...
                 (+a^2+b^2-c^2)/(2*a*b)];
        
        % If all cosines are non-negative, all angles are less than pi/2
        if 0 <= min(cos_T)
            % Find the circumradius
            RR(i) = a/(2*sqrt(1-cos_T(1)^2));
        end
        % If there is a negative cosine, we don't register the circumradius
    end
    % Find the maximum
    i_max = find(RR==max(RR));
    T = TT(i_max(1),:);
end