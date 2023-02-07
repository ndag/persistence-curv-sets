%files = ["D41_S1", "D41_S2", "D41_Torus"];
%dates = ["02_27_2020", "02_27_2020", "02_28_2020"];
%max_pers = [1, 1, 1];

%files = ["D41_S1"];
%dates = ["02_27_2020"];
%max_pers = [2];

files = ["D41_Torus"];
dates = ["02_28_2020"];
max_pers = [1];

norms = [1, 1, 1];
save_to_file = false;

for i=1:numel(files)
    % Load results
    file_name = sprintf("%s_%s", files(i), dates(i));
    load(sprintf("%s.mat", file_name));
    
    % If we are normalizing, we load the norm
    norm = norms(i);
    
%     % Set the color scheme for the graph
%     % We color points according to their distance to the diagonal
%     colorMap = jet(length(unique(z)));
%     set(gcf, 'ColorMap', colorMap);
%     h = plot(x,y);
    
    % Put the results in one matrix
    bd_times = cell2mat(bd_results);
    xx = bd_times(:,1)/norm;
    yy = bd_times(:,2)/norm;

    % Calculate the persistence of each point and get colors to paint
    % according to it
    pers = (yy-xx)/2;
    %rgb = vals2colormap(pers, 'jet', [0, max_pers(i)/norm]);
    
    % Plot and color according to the persistence of each point
    figure();
    colormap('jet');
    scatter(xx,yy, 1.5, pers, 'Filled');
    %colorbar;
    xlim([0, max(xx)]);
    ylim([0, max(yy)]);
    
    if save_to_file
        file_name = sprintf('%s/%s.png', pwd, files(i));
        saveas(gcf, file_name);
    end
end
