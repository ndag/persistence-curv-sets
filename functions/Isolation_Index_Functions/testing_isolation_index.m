%% Test 1:
% A square in S^1
X = [0, 0.5, 1, 1.5];
dm = abs(X-X');
dm = min(dm, 2-dm);

% I expect the splits:
% - [1,2], [3,4]
% - [2,3], [4,1],
% each with isolation index 1.
[alpha, splits] = IsolationIndexMatrix(dm);

N_splits = size(splits,1);
for i=1:N_splits
    fprintf('----------------\n')
    disp(alpha(i));
    disp(splits{i,1});
    disp(splits{i,2});
end

%% Test 2:
% A tree with 6 points
A = [0, 0, 1, 0, 0, 0;
     0, 0, 1, 0, 0, 0;
     1, 1, 0, 1, 0, 0;
     0, 0, 1, 0, 1, 1;
     0, 0, 0, 1, 0, 0;
     0, 0, 0, 1, 0, 0];

close all;
plot(graph(A));

dm = [0, 2, 1, 2, 3, 3;
      2, 0, 1, 2, 3, 3;
      1, 1, 0, 1, 2, 2;
      2, 2, 1, 0, 1, 1;
      3, 3, 2, 1, 0, 2;
      3, 3, 2, 1, 2, 0];

% I expect the splits:
% - [1], setdiff(1:6,1)
% - [2], setdiff(1:6,2)
% - [5], setdiff(1:6,5)
% - [6], setdiff(1:6,6)
% - [1,2,3], [4,5,6]
% All with isolation index 1
tic
[alpha, splits] = IsolationIndexMatrix(dm);
toc

N_splits = size(splits,1);
for i=1:N_splits
    fprintf('----------------\n')
    fprintf('alpha = %d\n', alpha(i));
    disp(splits{i,1});
    disp(splits{i,2});
end