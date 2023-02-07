clear; close all; clc;

addpath('../')
addpath('./export_fig')
tic

idx_0 = 2;
results_folder_0 = sprintf('../results/Curves_Set_%i', idx_0);

save_image = true;

% Load results
file_name_0 = sprintf('%s/Curves_%i', results_folder_0, idx_0);

load(sprintf('%s.mat', file_name_0));

%     % Set the color scheme for the graph
%     % We color points according to their distance to the diagonal
%     colorMap = jet(length(unique(z)));
%     set(gcf, 'ColorMap', colorMap);
%     h = plot(x,y);

%% Plot the results
for c=1:nCurves
% for c=1:1
    VV = Vertices{c};
    PP = Curves{c};
    
    PDs = BD_Curves{c};
    bd_times = unique(BD_times{c}, 'rows');
    
    % See the results
    if save_image
        f = figure('visible', 'off');
    else
        f = figure('visible', 'on');
    end
    
    % Title
    name = sprintf('Curve %i', c);
    sgtitle(name);
    
    % --------------------------------------
    % Plot the curve and the sampled points
    % --------------------------------------
    subplot(1,3,1);
    hold on;
    plot(VV(:,1), VV(:,2));
    plot(PP(:,1), PP(:,2),'r.');
    hold off;
    
    % --------------------------------------
    % Plot the persistence diagram
    % --------------------------------------
    subplot(1,3,2);
    hold on;
    plotpersistencediagrams(PDs);
    
    % --------------------------------------
    % Plot the persistence set
    % --------------------------------------
    subplot(1,3,3);
    %%
    hold on;
    
    % Main diagonal
    tt=0:0.01:max(bd_times(:));
    plot(tt,tt,'red');
    
    colormap('jet');
    scatplot(bd_times(:,1),bd_times(:,2),'circles');

    % Remap the colorbar to reflect the density 
    c = colorbar;
    % c.TickLabels = num2cell(linspace(1,2,6))
    % colorbar('delete')
    
    xlim([0, max(bd_times(:,1))]);
    ylim([0, max(bd_times(:,2))]);

    % set(gca, 'fontsize', 16)
    % set(gca, 'XTick', (0:0.5:3));
    
    hold off;
    
    % --------------------------------------
    % Save the results to a file
    % --------------------------------------
    if save_image
        % Reshape the canvas
        x_width=12;  y_width=4;
        set(gcf, 'PaperPosition', [0 0 x_width y_width]);
        
        saveas(gcf, sprintf('%s/%s.png', results_folder_0, name));  %#ok<UNRCH>
        savefig(sprintf('%s/%s.fig', results_folder_0, name));  %#ok<UNRCH>
    end
end % Of one of the nCurves repetitions

% if save_image
%     file_name = sprintf('%s_measure.pdf', files(i));
%     set(gcf, 'Color', 'w')
%     export_fig(file_name, '-q0')
%     
%     print(gcf,'-depsc', '-painters',file_name);
%     % Clean the output file with an external tool
%     epsclean(file_name);
% end
toc