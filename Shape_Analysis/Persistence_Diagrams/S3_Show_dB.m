% Show dB matrices
close all

% Reload matrices
dBs = zeros(maxHomDim+1, nchoosek(nShapes,2));
for d=1:maxHomDim+1
    mFile2 = matfile(sprintf('%s/Distances/dB_%i.mat', results_folder,d-1));
    dBs(d,:) = mFile2.dB;
end
mFile2 = matfile(sprintf('%s/Distances/dB_max.mat', results_folder), 'Writable', true);
C = mFile2.dB;

for d=1:maxHomDim+1
    f = figure('visible', 'off');
    colormap('hot');
    imagesc(squareform(dBs(d,:)));
    hold on
    title(sprintf('dB_%i', d-1), 'interpreter', 'none')
    colorbar;
    
    axis('equal')
    xlim([1,nShapes]);
    ylim([1,nShapes]);
    
    set(gca, 'xtick', labCenter, 'xticklabel', labels)
    set(gca, 'ytick', labCenter, 'yticklabel', labels)
    
    set(f, 'Color', 'w')
    export_fig(f, sprintf('%s/Distances/dB_%i.pdf', results_folder, d-1))
end

f = figure('visible', 'off');
colormap('hot');
imagesc(squareform(C));
hold on
title('dB_max', 'interpreter', 'none')
colorbar;

axis('equal')
xlim([1,nShapes]);
ylim([1,nShapes]);

set(gca, 'xtick', labCenter, 'xticklabel', labels)
set(gca, 'ytick', labCenter, 'yticklabel', labels)

set(f, 'Color', 'w')
export_fig(f, sprintf('%s/Distances/dB_max.pdf', results_folder));
