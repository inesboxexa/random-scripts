

mean127 = mean(data127,'omitnan'); mean127(isnan(mean127)) = []; 
mean145 = mean(data145,'omitnan'); mean145(isnan(mean145)) = [];
mean163 = mean(data163,'omitnan'); mean163(isnan(mean163)) = [];
mean181 = mean(data181,'omitnan'); mean181(isnan(mean181)) = [];

titles_array = [ "Glutamine [M-H3O]-", "Glutamine [M-H]-", "Glutamine [M+OH]-", "Glutamine [M+Cl]-" ]; 
    
figure
adduct = 0;
for data = { mean127, mean145, mean163, mean181 }
    adduct = adduct+1;
    
    wt = data{1}(:,1:13)';
    kras = data{1}(:,14:26)';
    apc = data{1}(:,27:39)';
    apckras = data{1}(:,40:53)';
    
    rows_num = max([size(wt,1),size(kras,1),size(apc,1),size(apckras,1)]);
    matrix = NaN*ones(rows_num,4);
    matrix(1:size(wt,1),1) = wt;
    matrix(1:size(kras,1),2) = kras;
    matrix(1:size(apc,1),3) = apc;
    matrix(1:size(apckras,1),4) = apckras;
    
    figure
    %subplot(2,4,adduct);
    
    for col = 1:size(matrix,2)
        title(titles_array(adduct))
        plot(-0.1+col,matrix(:,col),'.','markersize',10,'color',[.5 .5 .5]); hold on;
        plot(0.1+col,mean(matrix(:,col),'omitnan'),'x','color',[.9 0 .5],'markersize',10,'linewidth',2)
        plot(0.1+col,mean(matrix(:,col),'omitnan')+std(matrix(:,col),'omitnan'),'v','color',[.9 0 .5],'markersize',5,'linewidth',1)
        plot(0.1+col,mean(matrix(:,col),'omitnan')-std(matrix(:,col),'omitnan'),'^','color',[.9 0 .5],'markersize',5,'linewidth',1)
    end
    
    axis([0 5 0.9*min(matrix(:)) 1.1*max(matrix(:))])
    ylabel('Mean Intensity')
    xticks(1:4)
    xticklabels({'WT','KRAS','APC','APC-KRAS'})
    
end