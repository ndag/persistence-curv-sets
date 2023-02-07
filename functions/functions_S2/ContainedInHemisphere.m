% Checks if a set of points X is contained in a hemisphere of S^2. To do
% that, I take several unit vectors v normal to a plane and check if X is
% on one side of the plane. They are if all dot products X(i,:)*v have the
% same sign.
function flag = ContainedInHemisphere(X, nAngles)
    if nargin<2
        nAngles = 25;
    end
    
    % Create a set with several unit normal vectors
    Theta = 0:2*pi/nAngles:2*pi;
    Phi   = 0:(pi/2)/nAngles:pi/2;
    
    [TT, PP] = meshgrid(Theta,Phi);
    TT = reshape(TT, [], 1);
    PP = reshape(PP, [], 1);
    NN = [cos(TT).*sin(PP), sin(TT).*sin(PP), cos(PP)];
    
    % Do the dot product of all these vectors with our set of points
    % Each row contains the dot product of one normal vector n=NN(i,:) with
    % all the points in X.
    test = NN*X';
    
    % If there is a row with all non-positive or non-negative values, then
    % X is contained in a hemisphere
    non_positive_amount = sum(test>=0,2);
    non_negative_amount = sum(test<=0,2);
    
    if ismember(6,non_positive_amount) || ismember(6,non_negative_amount)
        flag=true;
    else
        flag=false;
    end
end
