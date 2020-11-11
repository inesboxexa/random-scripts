
%% Drug ion images

cd('X:\Alex\MS\2020090400_NiCEACQ_AD_DT_JZ_NB_TM\Data\SIMS data\Analysed\Sample 1_2\')
load('drugImagesList.mat')

%

original_folder = 'X:\Alex\MS\2020090400_NiCEACQ_AD_DT_JZ_NB_TM\Data\MVA\Sample 1_2\Target region\';
cd(original_folder)

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

%% Combining drugs

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

%%

original_folder = 'X:\Alex\MS\2020090400_NiCEACQ_AD_DT_JZ_NB_TM\Data\MVA\Sample 1_2\Target region\';

cd(original_folder)
load('combined_drugs_binned.mat')

bottom = combined_drugs_binned;

%

original_folder = 'X:\Alex\MS\2020090400_NiCEACQ_AD_DT_JZ_NB_TM\Data\MVA\Sample 1_1\';

cd(original_folder)
load('combined_drugs_binned.mat')

top = combined_drugs_binned;

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

%%

titles = ["1670 2nd Harmonic", "1670 SRS", "2850 CARS", "2850 2nd Harmonic", "2850 SRS",  "TPEF 2nd Harmonic" ];

%

red_th = 0.7;
green_th = 0.7;
blue_th = 0.2;

original_folder = 'X:\Alex\MS\2020090400_NiCEACQ_AD_DT_JZ_NB_TM\Data\MVA\Sample 1_2\Target region\';

cd(original_folder)
load('combined_drugs_binned.mat')
load('RAMANimage')
load('raman_skin_mask_kmean3')

skin_mask = reshape(raman_skin_mask_kmean3,size(RAMANimage,1),size(RAMANimage,2));
redi = 3;   
red_bottom = RAMANimage(:,:,redi).*skin_mask;
greeni = 4; 
green_bottom = RAMANimage(:,:,greeni).*skin_mask;
blue_bottom = combined_drugs_binned.*skin_mask;

%

original_folder = 'X:\Alex\MS\2020090400_NiCEACQ_AD_DT_JZ_NB_TM\Data\MVA\Sample 1_1\';

cd(original_folder)
load('combined_drugs_binned.mat')
load('RAMANimage')
load('raman_skin_mask_kmean3')

skin_mask = reshape(raman_skin_mask_kmean3,size(RAMANimage,1),size(RAMANimage,2));
redi = 3;   
red_top = RAMANimage(:,:,redi).*skin_mask;
greeni = 4; 
green_top = RAMANimage(:,:,greeni).*skin_mask;
blue_top = combined_drugs_binned.*skin_mask;

%

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

