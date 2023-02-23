% clear; close all; clc;

dHs = zeros(numel(Ks), nchoosek(nShapes,2));
dWs = zeros(numel(Ks), nchoosek(nShapes,2));
for idx=1:numel(Ks)
    k=Ks(idx);
    n=2*k+2;
    
    if useShortlist
        dH_file = matfile(sprintf('%s/Distances/dH_%i_%i_(%i).mat', results_folder,n,k,numel(short)));
        dW_file = matfile(sprintf('%s/Distances/dW_%i_%i_(%i).mat', results_folder,n,k,numel(short)));
    else
        dH_file = matfile(sprintf('%s/Distances/dH_%i_%i.mat', results_folder,n,k));
        dW_file = matfile(sprintf('%s/Distances/dW_%i_%i.mat', results_folder,n,k));
    end
    dH = dH_file.dH;
    dW = dW_file.dW;
    
    if size(dH,1)>nchoosek(numel(names),2)
        dH = squareform(dH);
        dW = squareform(dW);
        dHs(idx,:) = squareform(dH(short,short));
        dWs(idx,:) = squareform(dW(short,short));
    else
        dHs(idx,:) = dH;
        dWs(idx,:) = dW;
    end
end

% Compute the maximum over k
dH_max = max(dHs,[],1);

% For dW, I need to weight each row by a binomial coefficient
% (it might be possible to improve that bound)
for ik=1:numel(Ks)
    k=Ks(ik);
    n=2*k+2;
    dWs(ik,:) = dWs(ik,:)/nchoosek(n,2);
end
dW_max = max(dWs,[],1);

% Save matrices
if useShortlist
    dH_file = matfile(sprintf('%s/Distances/dH_max_(%i).mat', results_folder,nShapes), 'Writable', true);
    dW_file = matfile(sprintf('%s/Distances/dW_max_(%i).mat', results_folder,nShapes), 'Writable', true);
else
    dH_file = matfile(sprintf('%s/Distances/dH_max.mat', results_folder), 'Writable', true);
    dW_file = matfile(sprintf('%s/Distances/dW_max.mat', results_folder), 'Writable', true);
end
dH_file.dH = dH_max;
dW_file.dW = dW_max;

%% Show heatmaps of the maximum Hausdorff matrix
f = figure('visible', 'off');

colormap('hot');
imagesc(squareform(dH_max));
hold on
title('dHausdorff -- max')
colorbar;

axis('equal')
xlim([1,nShapes]);
ylim([1,nShapes]);

set(gca, 'xtick', labCenter, 'xticklabel', labels)
set(gca, 'ytick', labCenter, 'yticklabel', labels)

% Stuff to remove whitespace
set(f, 'PaperUnits','centimeters');
set(f, 'Units','centimeters');
pos=get(f,'Position');
set(f, 'PaperSize', [pos(3) pos(4)]);
set(f, 'PaperPositionMode', 'manual');
set(f, 'PaperPosition',[0 0 pos(3) pos(4)]);
if useShortlist
    saveas(f, sprintf('%s/Distances/dH_max_(%i).pdf', results_folder,nShapes));
else
    saveas(f, sprintf('%s/Distances/dH_max.pdf', results_folder));
end

%% Show heatmaps of the maximum Hausdorff matrix
f = figure('visible', 'off');

colormap('hot');
imagesc(squareform(dW_max));
hold on
title('dWasserstein -- max')
colorbar;

axis('equal')
xlim([1,nShapes]);
ylim([1,nShapes]);

set(gca, 'xtick', labCenter, 'xticklabel', labels)
set(gca, 'ytick', labCenter, 'yticklabel', labels)

% Stuff to remove whitespace
set(f, 'Color', 'w')
if useShortlist
    export_fig(f, sprintf('%s/Distances/dW_max_(%i).pdf', results_folder,nShapes));
else
    export_fig(f, sprintf('%s/Distances/dW_max.pdf', results_folder));
end
