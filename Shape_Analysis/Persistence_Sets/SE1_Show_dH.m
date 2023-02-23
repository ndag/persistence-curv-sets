% clear; close all; clc;

dHs = zeros(numel(Ks), nchoosek(nShapes,2));
for idx=1:numel(Ks)
    k=Ks(idx);
    n=2*k+2;
    
    if useShortlist
        dH_file = matfile(sprintf('%s/Distances/dH_%i_%i_(%i).mat', results_folder,n,k,numel(short)));
    else
        dH_file = matfile(sprintf('%s/Distances/dH_%i_%i.mat', results_folder,n,k));
    end
    dH = dH_file.dH;
    
    if size(dH,1)>nchoosek(nShapes,2)
        dH = squareform(dH);
        dHs(idx,:) = squareform(dH(short,short));
    else
        dHs(idx,:) = dH;
    end
end

%% Show heatmaps of the Hausdorff matrices
% p = numSubplots(numel(Ks));

for idx=1:numel(Ks)
    f = figure('visible', 'off');
    
    k=Ks(idx);
    n=2*k+2;
    dH=dHs(idx,:);
    
    % subplot(p(1),p(2), idx);
    colormap('hot');
    imagesc(squareform(dH));
    hold on
    title(sprintf('dHausdorff D_{%i,%i}',n,k))
    colorbar;
    
    axis('equal')
    xlim([1,nShapes]);
    ylim([1,nShapes]);
    
    set(gca, 'xtick', labCenter, 'xticklabel', labels)
    set(gca, 'ytick', labCenter, 'yticklabel', labels)
    
    % Stuff to remove whitespace
    set(f, 'Color', 'w')
    if useShortlist
        export_fig(f, sprintf('%s/Distances/dH_%i_%i_(%i).pdf', results_folder,n,k,nShapes));
    else
        export_fig(f, sprintf('%s/Distances/dH_%i_%i.pdf', results_folder,n,k));
    end
end
