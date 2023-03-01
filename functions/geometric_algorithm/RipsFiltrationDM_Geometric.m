%Programmer: Chris Tralie
%Purpose: A wrapper around Uli Bauer's ripser
% Modified by Mario Gomez
% Modifications:
%           - Changed the code to extract the upper diagonal of the matrix.
%           - We now write the matrix to a binary file instead of a txt.
%Parameters: D (N x N distance matrix)
%           maxHomDim: The maximum homological dimension (default 1)
%           thresh: The maximum length at which to add an edge
%           coeff: The field coefficients to use (default 2)

function [PDs] = RipsFiltrationDM_Mario(D, maxHomDim, thresh, coeff)
    if nargin < 3
        thresh = max(D(:))*2;
    end
    if nargin < 4
        coeff = 2;
    end
    
    % Write the matrix to a csv file so that Ripser can read it
    file_name = sprintf('temp/DLower.csv');
    writematrix(D, file_name);

    % Call Ripser throught the terminal
    [~, cmdout] = system(sprintf('top-dim-ph --dim %i --modulus %i --threshold %f --format distance %s', maxHomDim, coeff, thresh, file_name));
    
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
