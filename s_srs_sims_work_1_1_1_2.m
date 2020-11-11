%%

load('X:\Alex\MS\2020090400_NiCEACQ_AD_DT_JZ_NB_TM\Data\MVA\Sample 1_1\combined_drugs_binned.mat');
top = combined_drugs_binned;

load('X:\Alex\MS\2020090400_NiCEACQ_AD_DT_JZ_NB_TM\Data\MVA\Sample 1_2\Target region\combined_drugs_binned.mat');
bottom = combined_drugs_binned;

%

image2plot_top = top;
image2plot_top(sum(image2plot_top,2)==0,:)=[];
image2plot_top(:,sum(image2plot_top,1)==0)=[];

image2plot_bottom = bottom;
image2plot_bottom(sum(image2plot_bottom,2)==0,:)=[];
image2plot_bottom(:,sum(image2plot_bottom,1)==0)=[];

image2plot = [
    [image2plot_top,    NaN*ones(size(image2plot_top,1),max(0,size(image2plot_bottom,2)-size(image2plot_top,2)))]
    [image2plot_bottom, NaN*ones(size(image2plot_bottom,1),max(0,size(image2plot_top,2)-size(image2plot_bottom,2)))]
    ];
    
fig0 = figure('units','normalized','outerposition',[0 0 1.0 1.0]);

imagesc(image2plot); axis image; axis off; colormap('hot'); colorbar;
title('drug ions combined', 'fontsize',16)

figname_char = 'drug & no drug combined drugs images.fig'; savefig(fig0,figname_char,'compact')
tifname_char = 'drug & no drug combined drugs images.tif'; saveas(fig0,tifname_char)

close all
clear fig0

%% Clustering Raman data

load('X:\Alex\MS\2020090400_NiCEACQ_AD_DT_JZ_NB_TM\Data\MVA\Sample 1_1\RAMANimage.mat');
RAMANimage_1 = RAMANimage;
load('X:\Alex\MS\2020090400_NiCEACQ_AD_DT_JZ_NB_TM\Data\MVA\Sample 1_1\raman_skin_mask_kmean3.mat');
mask1 = raman_skin_mask_kmean3;
load('X:\Alex\MS\2020090400_NiCEACQ_AD_DT_JZ_NB_TM\Data\MVA\Sample 1_1\combined_drugs_binned.mat');
combined_drugs_binned1 = combined_drugs_binned;

load('X:\Alex\MS\2020090400_NiCEACQ_AD_DT_JZ_NB_TM\Data\MVA\Sample 1_2\Target region\RAMANimage.mat');
RAMANimage_2 = RAMANimage;
load('X:\Alex\MS\2020090400_NiCEACQ_AD_DT_JZ_NB_TM\Data\MVA\Sample 1_2\Target region\raman_skin_mask_kmean3.mat');
mask2 = raman_skin_mask_kmean3;
load('X:\Alex\MS\2020090400_NiCEACQ_AD_DT_JZ_NB_TM\Data\MVA\Sample 1_2\Target region\combined_drugs_binned.mat');
combined_drugs_binned2 = combined_drugs_binned;


data1 = reshape(RAMANimage_1,[],size(RAMANimage_1,3));
data2 = reshape(RAMANimage_2,[],size(RAMANimage_2,3));

for K = 2:6

    [ idx10, ~, ~, ~ ] = kmeans(data1(mask1,:), K, 'distance', 'correlation');
    idx1 = zeros(size(mask1,1),1);
    idx1(mask1,:) = idx10;
    
    [ idx20, ~, ~, ~ ] = kmeans(data2(mask2,:), K, 'distance', 'correlation');
    idx2 = zeros(size(mask2,1),1);
    idx2(mask2,:) = idx20;
    
    idx2(idx2>0) = idx2(idx2>0)+max(idx1);
    
    idx12D = reshape(idx1,size(RAMANimage_1,1),size(RAMANimage_1,2));
    idx22D = reshape(idx2,size(RAMANimage_2,1),size(RAMANimage_2,2));
    
    image2plot = [
        [idx12D,    NaN*ones(size(idx12D,1),max(0,size(idx22D,2)-size(idx12D,2)))]
        [idx22D,    NaN*ones(size(idx22D,1),max(0,size(idx12D,2)-size(idx22D,2)))]
        ];
    image2plot(sum(image2plot,2)==0,:)=[];
    image2plot(:,sum(image2plot,1)==0)=[];

    fig0 = figure('units','normalized','outerposition',[0 0 .7 .7]);
    cmap = [0 0 0; tab10(2*K)];
    imagesc(image2plot); axis image; axis off; colormap(cmap); colorbar
    title(['k-means clustering (clusters number: ' num2str(K) ' , distance metric: correlation)'],'fontsize',20)
    
    sub_folder = ['skin only raman data kmeans ', num2str(K), ' clusters'];
    
    mkdir(sub_folder)
    cd(sub_folder)
    
    save('idx1','idx1')
    save('idx2','idx2')
    
    figname_char = 'drug and no drug all clusters image.fig'; savefig(fig0,figname_char,'compact')
    tifname_char = 'drug and no drug all clusters image.tif'; saveas(fig0,tifname_char)
    
    idx = [
        idx1
        idx2
        ];
    
    combined_drugs_binned = [
        reshape(combined_drugs_binned1,[],1)
        reshape(combined_drugs_binned2,[],1)
        ];
    
    fig1 = figure('units','normalized','outerposition',[0 0 .3 .3]);
    for k=1:max(idx)
        data4box = combined_drugs_binned.*(idx==k);
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
    
    figname_char = 'drug and no drug bar plots combined drugs.fig'; savefig(fig1,figname_char,'compact')
    tifname_char = 'drug and no drug bar plots combined drugs.tif'; saveas(fig1,tifname_char)    
    
    close all
    clear fig0 fig1
    
    cd(original_folder)
    
end 


%% Overlay Drug and not drugged tissues

redi = 3; 
greeni = 4;

cd('X:\Alex\MS\2020090400_NiCEACQ_AD_DT_JZ_NB_TM\Data\MVA\Sample 1_1\')
load('combined_drugs_binned.mat')
load('RAMANimage')
load('raman_skin_mask_kmean3')

skin_mask = reshape(raman_skin_mask_kmean3,size(RAMANimage,1),size(RAMANimage,2));
red_top = RAMANimage(:,:,redi).*skin_mask;
green_top = RAMANimage(:,:,greeni).*skin_mask;
blue_top = combined_drugs_binned.*skin_mask;

cd('X:\Alex\MS\2020090400_NiCEACQ_AD_DT_JZ_NB_TM\Data\MVA\Sample 1_2\Target region\')
load('combined_drugs_binned.mat')
load('RAMANimage')
load('raman_skin_mask_kmean3')

skin_mask = reshape(raman_skin_mask_kmean3,size(RAMANimage,1),size(RAMANimage,2));
red_bottom = RAMANimage(:,:,redi).*skin_mask;
green_bottom = RAMANimage(:,:,greeni).*skin_mask;
blue_bottom = combined_drugs_binned.*skin_mask;

%

red_th = 0.8;
green_th = 0.8;
blue_th = 0.3;

red = [
    [red_top,    NaN*ones(max(0,size(red_bottom,2)-size(red_top,2)))]
    [red_bottom, NaN*ones(max(0,size(red_top,2)-size(red_bottom,2)))]
    ];

green = [
    [green_top,    NaN*ones(max(0,size(green_bottom,2)-size(green_top,2)))]
    [green_bottom, NaN*ones(max(0,size(green_top,2)-size(green_bottom,2)))]
    ];

blue = [
    [blue_top,    NaN*ones(max(0,size(blue_bottom,2)-size(blue_top,2)))]
    [blue_bottom, NaN*ones(max(0,size(blue_top,2)-size(blue_bottom,2)))]
    ];

%

RGB = red+green+blue;
red = (red-min(red(:)))./max(red(:)-min(red(:)));
red = imadjust(red,[0 red_th]);
red(sum(RGB,2)==0, :) = [];
red(:, sum(RGB,1)==0) = [];
green = (green-min(green(:)))./max(green(:)-min(green(:)));
green = imadjust(green,[0 green_th]);
green(sum(RGB,2)==0, :) = [];
green(:, sum(RGB,1)==0) = [];
blue = (blue-min(blue(:)))./max(blue(:)-min(blue(:)));
blue = imadjust(blue,[0 blue_th]);
blue(sum(RGB,2)==0, :) = [];
blue(:, sum(RGB,1)==0) = [];

image2plot = NaN*ones(size(red,1),size(red,2),3);
image2plot(:,:,1) = red;
image2plot(:,:,2) = green;
image2plot(:,:,3) = blue;

fig0 = figure('units','normalized','outerposition',[0 0 1.0 1.0]);
imshow(image2plot, 'InitialMagnification', 400); axis image; axis off;
title({[char(titles(redi)), ' (red) - ', char(titles(greeni)), ' (green) - drug (blue)']},'fontsize',18)

figname_char = ['drug & no drug overlay ' char(titles(redi)) ' combined drug.fig']; savefig(fig0,figname_char,'compact')
tifname_char = ['drug & no drug overlay ' char(titles(redi)) ' combined drug.tif']; saveas(fig0,tifname_char)

close all
clear fig0

%%

redi = 5; 
greeni = 4;

cd('X:\Alex\MS\2020090400_NiCEACQ_AD_DT_JZ_NB_TM\Data\MVA\Sample 1_1\')
load('combined_drugs_binned.mat')
load('RAMANimage')
load('raman_skin_mask_kmean3')

skin_mask = reshape(raman_skin_mask_kmean3,size(RAMANimage,1),size(RAMANimage,2));
red_top = RAMANimage(:,:,redi).*skin_mask;
green_top = RAMANimage(:,:,greeni).*skin_mask;
blue_top = combined_drugs_binned.*skin_mask;

cd('X:\Alex\MS\2020090400_NiCEACQ_AD_DT_JZ_NB_TM\Data\MVA\Sample 1_2\Target region\')
load('combined_drugs_binned.mat')
load('RAMANimage')
load('raman_skin_mask_kmean3')

skin_mask = reshape(raman_skin_mask_kmean3,size(RAMANimage,1),size(RAMANimage,2));
red_bottom = RAMANimage(:,:,redi).*skin_mask;
green_bottom = RAMANimage(:,:,greeni).*skin_mask;
blue_bottom = combined_drugs_binned.*skin_mask;

%

red_th = 0.8;
green_th = 0.8;
blue_th = 0.3;

red = [
    [red_top,    NaN*ones(max(0,size(red_bottom,2)-size(red_top,2)))]
    [red_bottom, NaN*ones(max(0,size(red_top,2)-size(red_bottom,2)))]
    ];

green = [
    [green_top,    NaN*ones(max(0,size(green_bottom,2)-size(green_top,2)))]
    [green_bottom, NaN*ones(max(0,size(green_top,2)-size(green_bottom,2)))]
    ];

blue = [
    [blue_top,    NaN*ones(max(0,size(blue_bottom,2)-size(blue_top,2)))]
    [blue_bottom, NaN*ones(max(0,size(blue_top,2)-size(blue_bottom,2)))]
    ];

%

RGB = red+green+blue;
red = (red-min(red(:)))./max(red(:)-min(red(:)));
red = imadjust(red,[0 red_th]);
red(sum(RGB,2)==0, :) = [];
red(:, sum(RGB,1)==0) = [];
green = (green-min(green(:)))./max(green(:)-min(green(:)));
green = imadjust(green,[0 green_th]);
green(sum(RGB,2)==0, :) = [];
green(:, sum(RGB,1)==0) = [];
blue = (blue-min(blue(:)))./max(blue(:)-min(blue(:)));
blue = imadjust(blue,[0 blue_th]);
blue(sum(RGB,2)==0, :) = [];
blue(:, sum(RGB,1)==0) = [];

image2plot = NaN*ones(size(red,1),size(red,2),3);
image2plot(:,:,1) = red;
image2plot(:,:,2) = green;
image2plot(:,:,3) = blue;

fig0 = figure('units','normalized','outerposition',[0 0 1.0 1.0]);
imshow(image2plot, 'InitialMagnification', 400); axis image; axis off;
title({[char(titles(redi)), ' (red) - ', char(titles(greeni)), ' (green) - drug (blue)']},'fontsize',18)

figname_char = ['drug & no drug overlay ' char(titles(redi)) ' combined drug.fig']; savefig(fig0,figname_char,'compact')
tifname_char = ['drug & no drug overlay ' char(titles(redi)) ' combined drug.tif']; saveas(fig0,tifname_char)

close all
clear fig0