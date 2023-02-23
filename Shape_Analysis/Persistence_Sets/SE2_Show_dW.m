% clear; close all; clc;

dWs = zeros(numel(Ks), nchoosek(nShapes,2));
for idx=1:numel(Ks)
    k=Ks(idx);
    n=2*k+2;
    
    if useShortlist
        dW_file = matfile(sprintf('%s/Distances/dW_%i_%i_(%i).mat', results_folder,n,k,numel(short)));
    else
        dW_file = matfile(sprintf('%s/Distances/dW_%i_%i.mat', results_folder,n,k));
    end
    dW = dW_file.dW;
    
    if size(dW,1)>nchoosek(nShapes,2)
        dW = squareform(dW);
        dWs(idx,:) = squareform(dW(short,short));
    else
        dWs(idx,:) = dW;
    end
end

%% Show heatmaps of the Wasserstein matrices
% p = numSubplots(numel(Ks));

for idx=1:numel(Ks)
    f = figure('visible', 'off');
    
    k=Ks(idx);
    n=2*k+2;
    dW=dWs(idx,:);
    
    % subplot(p(1),p(2), idx);
    colormap('hot');
    imagesc(squareform(dW));
    hold on
    title(sprintf('dWasserstein U_{%i,%i}',n,k))
    colorbar;
    
    axis('equal')
    xlim([1,nShapes]);
    ylim([1,nShapes]);
    
    set(gca, 'xtick', labCenter, 'xticklabel', labels)
    set(gca, 'ytick', labCenter, 'yticklabel', labels)
    
    % Stuff to remove whitespace
    set(f, 'Color', 'w')
    if useShortlist
        export_fig(f, sprintf('%s/Distances/dW_%i_%i_(%i).pdf', results_folder,n,k,nShapes));
    else
        export_fig(f, sprintf('%s/Distances/dW_%i_%i.pdf', results_folder,n,k));
    end
end
