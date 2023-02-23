%Programmer: Chris Tralie
%Purpose: A wrapper around Uli Bauer's ripser.
% Modified by Mario Gomez
% Modifications: This function is an adaptation of RipsFiltrationDM that
% can be called inside a parallel loop. In order for Ripser to access the
% distance matrix, we have to write it to disk but, in a parallel
% environment, several workers can call this function at the same time.
% In order to avoid file conflicts, we assign a separate file to each
% worker.
%Parameters: D (N x N distance matrix)
%           maxHomDim: The maximum homological dimension (default 1)
%           thresh: The maximum length at which to add an edge
%           coeff: The field coefficients to use (default 2)
%           t: Index of the worker calling RipsFiltrationDMPar
%           (you can get it with t=GetCurrentTask();)
function [PDs] = RipsFiltrationDMPar(D, maxHomDim, t, name, thresh, coeff)
    if nargin < 4
        name="";
    end
    if nargin < 5
        thresh = max(D(:))*2;
    end
    if nargin < 6
        coeff = 2;
    end
    
    % Step 1: Extract and output lower triangular distance matrix
    n = length(D);
    DLT = D(boolean(triu(ones(n)-eye(n))));
    
    % Step 2: Write the matrix to a binary file so that Ripser can read it
    file_name = sprintf('temp/DLower_%s_%02i.bin',name,t);
    fileID = fopen(file_name,'w');
    fwrite(fileID, DLT, 'float', 'l');
    fclose(fileID);

    % Call Ripser throught the terminal
    [~, cmdout] = system(sprintf('~/MATLAB/functions/ripser/ripser-coeff --dim %i --threshold %g --modulus %i --format binary %s', maxHomDim, thresh, coeff, file_name));
    
    lines = strsplit(cmdout, '\n');
    PDs = {};
    for ii = 1:length(lines)
        if length(lines{ii}) < 4
            continue
        end
        if strcmp(lines{ii}(1:4), 'dist')
            continue
        elseif strcmp(lines{ii}(1:4), 'valu')
            continue
        elseif strcmp(lines{ii}(1:4), 'pers')
            PDs{end+1} = {};
        elseif strcmp(lines{ii}(1:4), 'spar')
            continue
        else
            s = strrep(lines{ii}, ' ', '');
            s = strrep(s, '[', '');
            s = strrep(s, ']', '');
            s = strrep(s, '(', '');
            s = strrep(s, ')', '');
            fields = strsplit(s, ',');
            PDs{end}{end+1} = str2double(fields);
        end
    end
    for ii = 1:length(PDs)
        PDs{ii} = cell2mat(PDs{ii}');
    end
end