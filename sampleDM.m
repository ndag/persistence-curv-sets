function [dm, I] = sampleDM(fn,n)
%load(fn,'dX')
X = randn(10,2);
nor = sum(X.^2,2).^.5;
X = X./(nor*[1 1]);
dX = L2_distance(X',X');
dX = 2*asin(dX/2);

I = [1 2 3 4];
dm = dX(I,I);
%I = randi(length(dX),n,1)               % get n indices
%dm = dX(I,I);
