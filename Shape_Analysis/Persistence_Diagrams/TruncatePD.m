function PD = TruncatePD(PD,t,keepInf)
if isempty(PD)
    return;
end

if nargin<3
    keepInf = true;
end

% Start by replacing all nans with infinity
PD(isnan(PD)) = inf;

% Replace all values larger than t with t
if keepInf
    PD(PD>t & isfinite(PD)) = t;
else
    PD(PD>t) = t;
end

% Remove the points that landed on the diagonal
PD( PD(:,1)==PD(:,2), :) = [];
end