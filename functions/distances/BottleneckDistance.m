% D1 and D2 are matrices of size N-by-2 and M-by-2 representing persistence
% diagrams
function dB = BottleneckDistance(D1, D2)
    % Write the diagrams to text files so that the binaries can read them
    file_1 = 'temp/D1.txt';
    file_2 = 'temp/D2.txt';
    writematrix(D1, file_1, 'Delimiter', ' ');
    writematrix(D2, file_2, 'Delimiter', ' ');

    % Call the binary code throught the terminal
    [~, cmdout] = system(sprintf('~/MATLAB/packages/bottleneck/bottleneck_dist %s %s',...
                                    file_1, file_2));
    
    dB = str2double(cmdout);
end