
K = 8;

clustmap = [ 0 0 0; tab10(K) ];

distributionsM0 = distributionsM;
smaller_masks_list0 = smaller_masks_list;

% Organize small masks labels

indicies4small_masks = [ 1 4 5 7 9 10 12 13 15 2 3 6 8 11 14 16 ];

distributionsM = distributionsM0(:,indicies4small_masks);
smaller_masks_list = smaller_masks_list0(indicies4small_masks,1);

% Scatter plots

% Percentages

kpercentage = zeros(length(unique(distributionsM(~isnan(distributionsM)))'),size(distributionsM,2));
for k = unique(distributionsM(~isnan(distributionsM)))'
    kpercentage(k,:) = sum(distributionsM == k,1)./sum(~isnan(distributionsM),1);    
end

x = repmat(1:size(kpercentage,2),size(kpercentage,1),1); x = x(:);
y = repmat((1:size(kpercentage,1))',1,size(kpercentage,2)); y = y(:);
sz = kpercentage; sz = sz(:).*5000; % sz = (sz-min(sz(sz>0))+1)./max(sz-min(sz(sz>0))+1);
c = repmat(clustmap(2:end,:),size(kpercentage,2),1);

% Plotting

figure;
hold on;
stem([0:size(kpercentage,2)+1]-0.5,size(kpercentage,2)*ones(1,size(kpercentage,2)+2),'w')
plot(0:size(kpercentage,2)+1,repmat([1:size(kpercentage,1)+1]'-0.5,1,size(0:size(kpercentage,2)+1,2)),'w')
scatter(x(sz>0),y(sz>0),sz(sz>0),c(sz>0,:),'Marker','.');
set(gca,'Color','k')
set(gca,'XColor','k','YColor','k','TickDir','out')
axis image;
ylim([0.5 size(kpercentage,1)+0.5]);
xticks(1:size(kpercentage,1))
xlim([0.5 size(kpercentage,2)+0.5]);
xticks(1:size(kpercentage,2))
for i = 1:size(smaller_masks_list,1); aux_xticklabels{i} = char(smaller_masks_list(i)); end
xticklabels(aux_xticklabels)
xtickangle(90)