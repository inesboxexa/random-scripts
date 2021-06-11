
% This script creates new small masks with background (to plot single ion
% images).

cd('X:\Beatson\Intracolonic tumour study\dpo\neg DESI\rois')

files = dir('*slide9*');

for i = 1:size(files,1)
    
    cd([ files(i).folder filesep files(i).name ])
    
    rois = dir('*-S9*');
    
    for ii = 1:size(rois,1)
       
        cd([ rois(ii).folder filesep rois(ii).name ])
        
        load('roi')
        
        figure; subplot(1,2,1); imagesc(roi.pixelSelection); axis image
        
        x1 = max(find(sum(roi.pixelSelection,1)>0,1,'first')-20,1);
        x2 = min(find(sum(roi.pixelSelection,1)>0,1,'last')+20,size(roi.pixelSelection,2));
        
        y1 = max(find(sum(roi.pixelSelection,2)>0,2,'first')-20,1);
        y2 = min(find(sum(roi.pixelSelection,2)>0,2,'last')+20,size(roi.pixelSelection,1));
        
        background_mask = zeros(roi.height,roi.width);
        background_mask(y1:y2,x1:x2) = 1;
                
        roi = RegionOfInterest(roi.width,roi.height);
        roi.addPixels(background_mask)
        
        subplot(1,2,2); imagesc(roi.pixelSelection); axis image
        
        cd('..')
        mkdir([ rois(ii).folder filesep rois(ii).name '-bg'])
        cd([ rois(ii).folder filesep rois(ii).name '-bg'])
        
        save('roi','roi')
        
    end
    
end