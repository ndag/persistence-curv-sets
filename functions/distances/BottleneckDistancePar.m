% D1 and D2 are matrices of size N-by-2 and M-by-2 representing persistence
% diagrams
% This Par version requires an extra argument t so that different workers
% don't read/write on the same file on disk
function dB = BottleneckDistancePar(D1, D2, t)
    % Write the diagrams to text files so that the binaries can read them
    file_1 = sprintf('temp/D1_%02i.txt',t);
    file_2 = sprintf('temp/D2_%02i.txt',t);
    writematrix(D1, file_1, 'Delimiter', ' ');
    writematrix(D2, file_2, 'Delimiter', ' ');

    % Call the binary code throught the terminal
    [~, cmdout] = system(sprintf('~/MATLAB/packages/bottleneck/bottleneck_dist %s %s',...
                                    file_1, file_2));
    
    dB = str2double(cmdout);
end