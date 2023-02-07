function adj_mats = read_graphs(file_name, use_text_file, save_to_file)
    if nargin==2
        save_to_file=false;
    end
    
    % Load the matrices from a text file
    if use_text_file
        % Read file and store each matrix into a cell entry (matrices should
        % be separated by an empty line in the text file)
        txt_from_file = fileread(sprintf('./data/%s.txt',file_name));
        adj_txt = strsplit(txt_from_file, '\n\n');
        
        adj_mats = cell(size(adj_txt));
        for i=1:numel(adj_txt)
            adj_mats{i} = str2num(adj_txt{i});
        end
    
    % Or load them directly from a .mat file
    else
        load(sprintf('./data/%s.mat', file_name), 'adj_mats');
    end
    
    if save_to_file
        save(sprintf('%s.mat',file_name),'adj_mats');
    end
end