
data_name = "pos api MALDI";

dataset_name = "positive AP MALDI tumour models";
main_mask = "tissue only";
mva_code_1 = "tsne 5 components";
norm_type = "no norm";
mva_code_2 = "4000 highest peaks";


regionsNumE = 1; % number of regions to keep
regionsNumI = 0; % number of regions to fill

imagei = 1;

file_name = "2020_01_28_50um_agcOff_Pos_tissue_image";

%         "2020_01_28_50um_agcOff_Pos_tissue_image"
%         "2020_01_28_50um_agcOff_Pos_tissueCDBF_image"
%         "2020_01_29_50um_agcOff_Pos_2nd_tissue_marcel_image"
%         "2020_01_29_50um_agcOff_Pos_tissue_marcel_image"

%

new_mask_name = "APC-normal"; % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% imagei = imagei+1;

roi_mask0 = 0;
roii = 1;

for component = [ 2 3 4 ] % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    load(['X:\Beatson\Intracolonic tumour study\dpo\' char(data_name) '\rois\' char(file_name) '\mva based rois\' char(dataset_name) '\' char(main_mask) '\' char(mva_code_1) '\' char(norm_type) '\mva ' char(mva_code_2) '\component ' num2str(component) '\roi.mat'])
    roi_mask0 = roi_mask0 + roi.pixelSelection;
    if roii; width = roi.width; height = roi.height; roii = 0; end
    clear roi
    
end

if regionsNumE > 0
    roi_mask1 = 0*roi_mask0;
else
    roi_mask1 = 1;
end

for regionsi = 1:regionsNumE
    figure(1000);
    roi_mask2 = roipoly(roi_mask0);
    roi_mask1 = roi_mask1 + roi_mask2;
    title({'Choose Areas to Keep.'})
end

roi_mask3 = 0*roi_mask0;
for regionsi = 1:regionsNumI
    figure(2000);
    roi_mask4 = roipoly(roi_mask0.*roi_mask1);
    roi_mask3 = roi_mask3 + roi_mask4;
    title({'Choose Areas to Fill In.'})
end

roi_mask = logical(roi_mask0 .* roi_mask1 + roi_mask3);

figure(3000);
imagesc(roi_mask); axis image; colormap gray;
title({'Final ROI.'})

roi = RegionOfInterest(width,height);
roi.addPixels(roi_mask)

mkdir(['X:\Beatson\Intracolonic tumour study\dpo\' char(data_name) '\rois\' char(file_name) '\' char(new_mask_name) '\'])
cd(['X:\Beatson\Intracolonic tumour study\dpo\' char(data_name) '\rois\' char(file_name) '\' char(new_mask_name) '\'])

save('roi','roi')

%

