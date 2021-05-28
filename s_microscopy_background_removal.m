
% This script reads the tif files in imagesPath, and saved new tif files
% with the background pixels set to white. The following steps are taken:
% (1) smothing of the image using a gaussian kernel
% (1 notes) size of the kernel can be adjusted - smaller kernels result in
% sharper clustering maps and in wholes in the tissue
% (2) k-means clustering using the smooth image
% (2 notes) number of clusters can be adjusted
% (3) selection of the background cluster
% (3 notes) at the moment, the largest clusters is set to be the
% background - you might have to change this depending on the image
% (4) increase the sharpness of the image
% (4 notes) at the moment, no parameters are adjusted
% (5) enhance the contrast
% (5 notes) at the moment, no parameters are adjusted
% (6) save the image with no background and enhanced contrast
% (6 notes) you could decide to save the sharper image instead

% By Teresa Murta, 28 May 2021


imagesPath = 'X:\MS\NICE ACR\H&E\H&E\10X\';
images = dir([imagesPath filesep '*.tif*']);

for filei = 1:size(images,1) % iterates through files
    
    t = Tiff([images(filei).folder filesep images(filei).name],'r');
    imageData = read(t); % reads image data (RGB)
    
    figure(filei)
    subplot(1,6,1); imshow(imageData); axis image; title('Original data')
    
    imageData_filter = imgaussfilt(imageData,8); % gaussian filter - change the number to change the size of the filter (smaller filters results in more wholes in the tissue and sharper edges)
    
    subplot(1,6,2); imshow(imageData_filter); axis image; title('Filtered data')
    
    data_array = double(reshape(imageData_filter,[],3));
    
    idx = kmeans(data_array,5,'Distance','cosine'); % k-means Clustering with 5 clusters
    
    subplot(1,6,3); imagesc(reshape(idx,size(imageData,1),size(imageData,2))); axis image off; title('Clustering map')
    
    k_array = [];
    for k = unique(idx)'
        k_array = [ k_array sum(idx==k) ];
    end
    
    [~, background_k] = max(k_array); % select the cluster with more pixels
    
    skinMask = uint8(reshape(idx~=background_k,size(imageData,1),size(imageData,2)));
    
    noBackgroundData = imageData;
    noBackgroundData(:,:,1) = skinMask.*imageData(:,:,1)+uint8(255*~skinMask);
    noBackgroundData(:,:,2) = skinMask.*imageData(:,:,2)+uint8(255*~skinMask);
    noBackgroundData(:,:,3) = skinMask.*imageData(:,:,3)+uint8(255*~skinMask);
        
    subplot(1,6,4); imshow(noBackgroundData); axis image; title('Background removed')
    
    noBackgroundData_sharp = imsharpen(noBackgroundData); % increases sharpeness
    
    subplot(1,6,5); imagesc(noBackgroundData_sharp); axis image off; title('Sharpeness enhanced')
    
    RGB = noBackgroundData;
    LAB = rgb2lab(RGB);

    L = LAB(:,:,1)/100;
    L = adapthisteq(L,'NumTiles',[8 8],'ClipLimit',0.005); % improves contrast
    LAB(:,:,1) = L*100;

    noBackgroundData_enhancedContrast = lab2rgb(LAB);
    
    subplot(1,6,6); imagesc(noBackgroundData_enhancedContrast); axis image off; title('Contrast enhanced')
    
    imwrite(noBackgroundData_enhancedContrast,[images(filei).name(1:end-4), '_nobk.tif'],'compression','none') % saves new tif file

end

%% Other tests - ignore

%     % Find the countour of the skin
%     
%     skinMask1 = 0*skinMask;
%     for rowi = 1:size(noBackgroundData,1)
%         aux = sum(noBackgroundData(rowi,:,1:2),3)~=255*2;
%         if sum(aux)>0
%             startX = find(sum(noBackgroundData(rowi,:,1:2),3)~=255*2,1,'first');
%             endX = find(sum(noBackgroundData(rowi,:,1:2),3)~=255*2,1,'last');
%             
%             skinMask1(rowi,startX:endX) = 1;
%         end
%     end
%     
%     skinMask2 = 0*skinMask;
%     for coli = 1:size(noBackgroundData,2)
%         aux = sum(noBackgroundData(:,coli,1:2),3)~=255*2;
%         if sum(aux)>0
%             startY = find(sum(noBackgroundData(:,coli,1:2),3)~=255*2,1,'first');
%             endY = find(sum(noBackgroundData(:,coli,1:2),3)~=255*2,1,'last');
%             
%             skinMask2(startY:endY,coli) = 1;
%         end
%     end
%     
%     skinMask = uint8(skinMask1.*skinMask2);
%             
%     noBackgroundData_special = imageData;
%     noBackgroundData_special(:,:,1) = skinMask.*imageData(:,:,1)+uint8(255*~skinMask);
%     noBackgroundData_special(:,:,2) = skinMask.*imageData(:,:,2)+uint8(255*~skinMask);
%     noBackgroundData_special(:,:,3) = skinMask.*imageData(:,:,3)+uint8(255*~skinMask);
%     
%     subplot(1,4,4); imshow(noBackgroundData_special); axis image
%     imwrite(noBackgroundData_special,[images(filei).name(1:end-4), '_bkremoved_special.tif'],'compression','none')



