function [] = plot_density(bd_times)
    colormap('jet');
    bd_times_live = bd_times( bd_times(:,1)>0, :);
    scatplot(bd_times_live(:,1),bd_times_live(:,2),'circles');
end
