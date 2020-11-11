
%%

original_folder = 'X:\Alex\MS\2020090400_NiCEACQ_AD_DT_JZ_NB_TM\Data\MVA\Sample 1_2\Target region\';

%%% Building 3D images

cd(original_folder)
load('combinedRegisteredData.mat')

RAMANimage = repmat(zeros(size(combinedMask)),1,1,size([ 2:6 8 ],2));
i = 0;
for zi = [ 2:6 8 ]
    i = i+1;
    z = zeros(size(combinedMask));
    z(logical(combinedMask)) = registeredRamanMSIonly(:,zi);
    z(isnan(z)) = 0;
    RAMANimage(:,:,i) = z;
end

save('RAMANimage','RAMANimage')

SIMSimage = repmat(zeros(size(combinedMask)),1,1,size(registeredMSIramanOnly,2));
for zi = 1:size(registeredMSIramanOnly,2)
    z = zeros(size(combinedMask));
    z(logical(combinedMask)) = registeredMSIramanOnly(:,zi);
    z(isnan(z)) = 0;
    SIMSimage(:,:,zi) = z;
end

%%% Drugs in SIMS data

load('X:\Alex\MS\2020090400_NiCEACQ_AD_DT_JZ_NB_TM\Data\SIMS data\Analysed\Sample 1_2\drugImagesList.mat')

%

% Plot drugs

fig0 = figure('units','normalized','outerposition',[0 0 1.0 1.0]);
SIMSdrugs = zeros(size(imagesList{1},1),size(imagesList{1},2),length(imagesList));
SIMSdrugs_binned = zeros(size(imagesList{1},1),size(imagesList{1},2),length(imagesList));
for i = 1:length(imagesList)
    
    image = imagesList{i};
    
    image2plot = image;
    image2plot(sum(image2plot,2)==0,:)=[];
    image2plot(:,sum(image2plot,1)==0)=[];
    
    subplot(4,3,i); imagesc(image2plot); axis image; axis off; colormap('hot'); colorbar;
    title(['m/z ', num2str(round(massToCreate(i),2))],'fontsize',16)
    
    SIMSdrugs(:,:,i) = imagesList{i};
    
    % Binning drug image
    
    square_size = 4;
    drug_image_binned = 0*image;
    for rowi = 1:(size(image,1)-square_size)
        for coli = 1:(size(image,2)-square_size)
            drug_image_binned(rowi+2,coli+2) = sum(sum(image(rowi:rowi+square_size,coli:coli+square_size)));
        end
    end
    
    SIMSdrugs_binned(:,:,i) = drug_image_binned;
    
    image2plot = drug_image_binned;
    image2plot(sum(image2plot,2)==0,:)=[];
    image2plot(:,sum(image2plot,1)==0)=[];
    
    subplot(4,3,i+6); imagesc(image2plot); axis image; axis off; colormap('hot'); colorbar;
    title(['m/z ', num2str(round(massToCreate(i),2)), ' (binned)'],'fontsize',16)
    
    image = image(:);
    image(image==0)=[];
    image(isnan(image))=[];
    peak_mean_counts(i) = mean(image);
end

figname_char = 'drug ions images.fig'; savefig(fig0,figname_char,'compact')
tifname_char = 'drug ions images.tif'; saveas(fig0,tifname_char)

close all
clear fig0

save('SIMSdrugs','SIMSdrugs')
save('SIMSdrugs_binned','SIMSdrugs_binned')

%%% Combining drugs

% Combined drug image

combined_drugs = sum(SIMSdrugs,3);

% Binning combined drug image

square_size = 4;
combined_drugs_binned = 0*combined_drugs;
for rowi = 1:(size(combined_drugs,1)-square_size)
    for coli = 1:(size(combined_drugs,2)-square_size)
        combined_drugs_binned(rowi+2,coli+2) = sum(sum(combined_drugs(rowi:rowi+square_size,coli:coli+square_size)));
    end
end

save('combined_drugs','combined_drugs')
save('combined_drugs_binned','combined_drugs_binned')

fig0 = figure('units','normalized','outerposition',[0 0 1.0 1.0]);

subplot(1,2,1)
image2plot = combined_drugs;
image2plot(sum(image2plot,2)==0,:)=[];
image2plot(:,sum(image2plot,1)==0)=[];
imagesc(image2plot); axis image; axis off; colormap('hot'); colorbar;
title('drug ions combined', 'fontsize',16)

subplot(1,2,2)
image2plot = combined_drugs_binned;
image2plot(sum(image2plot,2)==0,:)=[];
image2plot(:,sum(image2plot,1)==0)=[];
imagesc(image2plot); axis image; axis off; colormap('hot'); colorbar;
title('drug ions combined (binned)','fontsize',16)

figname_char = 'combined drugs images.fig'; savefig(fig0,figname_char,'compact')
tifname_char = 'combined drugs images.tif'; saveas(fig0,tifname_char)

close all
clear fig0

%%%

% Spatial segmentation of Raman data

dataM = RAMANimage;

data = reshape(dataM,[],size(dataM,3));
mask = sum(data,2)==0;
data(mask,:)=[];

K = 3;

[ idx0, ~, ~, ~ ] = kmeans(data,K,'distance','correlation');

idx = zeros(size(mask,1),1);
idx(~mask,:) = idx0;

fig0 = figure('units','normalized','outerposition',[0 0 .7 .7]);
image2plot = reshape(idx,size(dataM,1),size(dataM,2));
imagesc(image2plot); axis image; axis off; colormap([0 0 0; tab10(K)]); colorbar
title(['k-means clustering (clusters number: ' num2str(K) ' , distance metric: correlation)'],'fontsize',16)
figname_char = 'all clusters image 2 define raman skin mask.fig'; savefig(fig0,figname_char,'compact')
tifname_char = 'all clusters image 2 define raman skin mask.tif'; saveas(fig0,tifname_char)

close all
clear fig0

%% Save Raman skin mask

raman_skin_mask_kmean3 = logical((idx==1)+(idx==3));

save('raman_skin_mask_kmean3','raman_skin_mask_kmean3')

%%

% Load Raman skin mask

load('raman_skin_mask_kmean3')

% Spatial segmentation of skin only Raman data

dataM = RAMANimage;

data = reshape(dataM,[],size(dataM,3));
mask = ~raman_skin_mask_kmean3;
data(mask,:)=[];

for K = 2:6
    
    [ idx0, ~, ~, ~ ] = kmeans(data, K, 'distance', 'correlation');
    
    idx = zeros(size(mask,1),1);
    idx(~mask,:) = idx0;
    
    fig0 = figure('units','normalized','outerposition',[0 0 .7 .7]);
    cmap = [0 0 0; tab10(K)];
    image2plot = reshape(idx,size(dataM,1),size(dataM,2));
    image2plot(sum(image2plot,2)==0,:)=[];
    image2plot(:,sum(image2plot,1)==0)=[];
    imagesc(image2plot); axis image; axis off; colormap(cmap); colorbar
    title(['k-means clustering (clusters number: ' num2str(K) ' , distance metric: correlation)'],'fontsize',20)
    
    sub_folder = ['skin only raman data kmeans ', num2str(K), ' clusters'];
    
    mkdir(sub_folder)
    cd(sub_folder)
    
    save('idx','idx')
    
    figname_char = 'all clusters image.fig'; savefig(fig0,figname_char,'compact')
    tifname_char = 'all clusters image.tif'; saveas(fig0,tifname_char)
    
    fig1 = figure('units','normalized','outerposition',[0 0 .1 .3]);
    data4box0 = reshape(combined_drugs_binned,[],1);
    for k=1:K
        data4box = data4box0.*(idx==k);
        data4box(data4box==0) = NaN;
        grid on
        b = bar(k, mean(data4box,'omitnan'));
        b.FaceColor = cmap(k+1,:);
        hold on
        errhigh = std(data4box,'omitnan')./sqrt(sum(~isnan(data4box))); % Standard error
        errlow  = std(data4box,'omitnan')./sqrt(sum(~isnan(data4box)));
        er = errorbar(k,mean(data4box,'omitnan'),errlow,errhigh);
        er.Color = [0 0 0];
        er.LineStyle = 'none';
        set(gca,'xticklabel',{[]})
    end
    ylabel('Intensity / a.u.')
    xlabel('Cluster')
    hold off
    
    figname_char = 'bar plots combined drugs.fig'; savefig(fig1,figname_char,'compact')
    tifname_char = 'bar plots combined drugs.tif'; saveas(fig1,tifname_char)
    
    %     fig2 = figure('units','normalized','outerposition',[0 0 .7 .7]);
    %     for i = 1:size(SIMSdrugs_binned,3)
    %
    %         data4box0 = reshape(SIMSdrugs_binned(:,:,i),[],1);
    %
    %         data4box = NaN*ones(size(data4box0,1),length(unique(idx(idx>0))));
    %         for ki=unique(idx(idx>0))'
    %             data4box(:,ki) = data4box0.*(idx==ki);
    %         end
    %         data4box(data4box==0) = NaN;
    %
    %         subplot(2,3,i)
    %         % boxplot(data4box)
    %         % plot(mean(data4box,'omitnan'),'xk')
    %         for k=1:K
    %             grid on
    %             b = bar(k, mean(data4box(:,k),'omitnan'));
    %             b.FaceColor = cmap(k+1,:);
    %             hold on
    %             errhigh = std(data4box(:,k),'omitnan')./sqrt(sum(~isnan(data4box(:,k)))); % Standard error
    %             errlow  = std(data4box(:,k),'omitnan')./sqrt(sum(~isnan(data4box(:,k))));
    %             er = errorbar(k,mean(data4box(:,k),'omitnan'),errlow,errhigh);
    %             er.Color = [0 0 0];
    %             er.LineStyle = 'none';
    %             set(gca,'xticklabel',{[]})
    %         end
    %         title(['m/z ', num2str(spectralChannels(drugindicies(i)))])
    %         ylabel('Intensity / a.u.')
    %         xlabel('Cluster')
    %         hold off
    %     end
    %
    %     figname_char = 'bar plots.fig'; savefig(fig2,figname_char,'compact')
    %     tifname_char = 'bar plots.tif'; saveas(fig2,tifname_char)
    
    close all
    clear fig0 fig1
    
    cd(original_folder)
    
end

%%%

% Plotting Raman channels

titles = ["1670 2nd Harmonic", "1670 SRS", "2850 CARS", "2850 2nd Harmonic", "2850 SRS",  "TPEF 2nd Harmonic" ];

fig0 = figure('units','normalized','outerposition',[0 0 1.0 1.0]);
for ramani = 1:size(RAMANimage,3)
    subplot(2,3,ramani);
    image2plot = RAMANimage(:,:,ramani);
    image2plot(sum(image2plot,2)==0,:)=[];
    image2plot(:,sum(image2plot,1)==0)=[];
    imagesc(image2plot); axis image; axis off; colormap('hot'); colorbar
    title({titles(ramani)},'fontsize',16)
end

figname_char = 'raman channels.fig'; savefig(fig0,figname_char,'compact')
tifname_char = 'raman channels.tif'; saveas(fig0,tifname_char)

close all
clear fig0

%%% Overlays

skin_mask = reshape(raman_skin_mask_kmean3,size(dataM,1),size(dataM,2));

%%% Overlays

skin_mask = reshape(raman_skin_mask_kmean3,size(dataM,1),size(dataM,2));

redi = 3;
red = RAMANimage(:,:,redi).*skin_mask;
greeni = 4;
green = RAMANimage(:,:,greeni).*skin_mask;

blue = combined_drugs_binned.*skin_mask;

RGB = red+green+blue;
red = (red-min(red(:)))./max(red(:)-min(red(:)));
red = imadjust(red,[0 0.7]);
red(sum(RGB,2)==0, :) = [];
red(:, sum(RGB,1)==0) = [];
green = (green-min(green(:)))./max(green(:)-min(green(:)));
green = imadjust(green,[0 .7]);
green(sum(RGB,2)==0, :) = [];
green(:, sum(RGB,1)==0) = [];
blue = (blue-min(blue(:)))./max(blue(:)-min(blue(:)));
blue = imadjust(blue,[0 0.2]);
blue(sum(RGB,2)==0, :) = [];
blue(:, sum(RGB,1)==0) = [];

image2plot = NaN*ones(size(red,1),size(red,2),3);
image2plot(:,:,1) = red;
image2plot(:,:,2) = green;
image2plot(:,:,3) = blue;

fig0 = figure('units','normalized','outerposition',[0 0 1.0 1.0]);
imshow(image2plot, 'InitialMagnification', 400); axis image; axis off;
title({[char(titles(redi)), ' (red) - ', char(titles(greeni)), ' (green) - drug (blue)']},'fontsize',18)
figname_char = ['overlay ' char(titles(redi)) ' combined drug.fig']; savefig(fig0,figname_char,'compact')
tifname_char = ['overlay ' char(titles(redi)) ' combined drug.tif']; saveas(fig0,tifname_char)

close all
clear fig0

redi = 5;
red = RAMANimage(:,:,redi).*skin_mask;
greeni = 4;
green = RAMANimage(:,:,greeni).*skin_mask;

blue = combined_drugs_binned.*skin_mask;

RGB = red+green+blue;
red = (red-min(red(:)))./max(red(:)-min(red(:)));
red = imadjust(red,[0 0.2]);
red(sum(RGB,2)==0, :) = [];
red(:, sum(RGB,1)==0) = [];
green = (green-min(green(:)))./max(green(:)-min(green(:)));
green = imadjust(green,[0 0.7]);
green(sum(RGB,2)==0, :) = [];
green(:, sum(RGB,1)==0) = [];
blue = (blue-min(blue(:)))./max(blue(:)-min(blue(:)));
blue = imadjust(blue,[0 0.2]);
blue(sum(RGB,2)==0, :) = [];
blue(:, sum(RGB,1)==0) = [];

image2plot = NaN*ones(size(red,1),size(red,2),3);
image2plot(:,:,1) = red;
image2plot(:,:,2) = green;
image2plot(:,:,3) = blue;

fig0 = figure('units','normalized','outerposition',[0 0 1.0 1.0]);
imshow(image2plot, 'InitialMagnification', 400); axis image; axis off;
title({[char(titles(redi)), ' (red) - ', char(titles(greeni)), ' (green) - drug (blue)']},'fontsize',18)
figname_char = ['overlay ' char(titles(redi)) ' combined drug.fig']; savefig(fig0,figname_char,'compact')
tifname_char = ['overlay ' char(titles(redi)) ' combined drug.tif']; saveas(fig0,tifname_char)

close all
clear fig0

% % Ion by ion
%
% for i = 1% 1:size(SIMSdrugs_binned,3)
%
%     ramani = 3;
%     red = RAMANimage(:,:,ramani).*skin_mask;
%     ramani = 4;
%     green = RAMANimage(:,:,ramani).*skin_mask;
%
%     blue = combined_drugs_binned.*skin_mask;
%     %blue = SIMSimage(:,:,drugindicies(i)).*skin_mask;
%
%     RGB = red+green+blue;
%     red = (red-min(red(:)))./max(red(:)-min(red(:)));
%     red = imadjust(red,[0 0.7]);
%     red(sum(RGB,2)==0, :) = [];
%     red(:, sum(RGB,1)==0) = [];
%     green = (green-min(green(:)))./max(green(:)-min(green(:)));
%     green = imadjust(green,[0 0.6]);
%     green(sum(RGB,2)==0, :) = [];
%     green(:, sum(RGB,1)==0) = [];
%     blue = (blue-min(blue(:)))./max(blue(:)-min(blue(:)));
%     blue = imadjust(blue,[0 0.3]);
%     blue(sum(RGB,2)==0, :) = [];
%     blue(:, sum(RGB,1)==0) = [];
%
%     image2plot = NaN*ones(size(red,1),size(red,2),3);
%     image2plot(:,:,1) = red;
%     image2plot(:,:,2) = green;
%     image2plot(:,:,3) = blue;
%
%     fig0 = figure('units','normalized','outerposition',[0 0 .7 .7]);
%     imshow(image2plot, 'InitialMagnification', 400); axis image; axis off;
%     title({[char(titles(6)), ' (red) - ', char(titles(2)), ' (green) - m/z ', char(num2str(spectralChannels(drugindicies(i)))), ' (blue)']},'fontsize',16)
%     figname_char = ['overlay ion ', char(num2str(spectralChannels(drugindicies(i)))) '.fig']; savefig(fig0,figname_char,'compact')
%     tifname_char = ['overlay ion ', char(num2str(spectralChannels(drugindicies(i)))) '.tif']; saveas(fig0,tifname_char)
%     %figname_char = ['overlay nonbinned ion ', char(num2str(spectralChannels(drugindicies(i)))) '.fig']; savefig(fig0,figname_char,'compact')
%     %tifname_char = ['overlay nonbinned ion ', char(num2str(spectralChannels(drugindicies(i)))) '.tif']; saveas(fig0,tifname_char)
%
%     close all
%     clear fig0
%
% end

