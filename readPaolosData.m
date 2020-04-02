folder = 'C:\Work\ICR_LCM\ICR\data\BR1458_5_cluster\';

cmz = h5read([folder 'count_matrix_multisample.h5'], '/cmz');
counts_mat = h5read([folder 'count_matrix_multisample.h5'], '/counts_mat');
dim_xy = h5read([folder 'count_matrix_multisample.h5'], '/dim_xy');

labels = h5read([folder 'mask.h5'], '/labels');
mask = h5read([folder 'mask.h5'], '/mask');

%%
counts = permute(reshape(counts_mat', dim_xy(1), dim_xy(2), []), [2 1 3]);

tissueMask = labels';
tissueMask = tissueMask == 1;

counts = reshape(counts, prod(dim_xy), []);
tissueMask = reshape(tissueMask, prod(dim_xy), []);

tissue = counts(tissueMask, :);

c = kmeans(tissue, 2, 'distance', 'cosine', 'replicates', 5);

clusterImage = zeros(size(tissueMask));
clusterImage(tissueMask) = c;

figure, imagesc(reshape(clusterImage, dim_xy(2), dim_xy(1)));

%%
tissueMask = clusterImage == 1;
tissueMask = reshape(tissueMask, prod(dim_xy), []);

tissue = counts(tissueMask, :);

c = kmeans(tissue, 5, 'distance', 'cosine', 'replicates', 5);

clusterImage = zeros(size(tissueMask));
clusterImage(tissueMask) = c;

k5 = reshape(clusterImage, dim_xy(2), dim_xy(1));

figure, imagesc(k5);

imwrite(ind2rgb(k5+1, parula(6)), 'C:\Work\ICR_LCM\ICR\data\BR1458_5_cluster\k=5.tif');

%%

minFeatureSize = 300;

kmeansImage = removeSmallClusters(k5, minFeatureSize);

k5rgb = ind2rgb(kmeansImage+1, parula(6));
figure, imagesc(k5rgb);
axis image; axis off;

imwrite(ind2rgb(kmeansImage+1, parula(6)), 'C:\Work\ICR_LCM\ICR\data\BR1458_5_cluster\k=5_removeSmallClusters.tif');

%%

newkmeansImage = zeros(size(kmeansImage));

for i = 1:max(kmeansImage(:))
    filled = imfill(kmeansImage == i, 'holes');
    
    newkmeansImage = newkmeansImage + double(newkmeansImage == 0) .* filled * i;
end

figure, imagesc(newkmeansImage);

roiList = RegionOfInterestList();
for i = 1:max(newkmeansImage(:))
    roi = RegionOfInterest(size(newkmeansImage, 2), size(newkmeansImage, 1));
    roi.addPixels(newkmeansImage == i);
    roi.setName(['k = ' num2str(i)]);
    roiList.add(roi);
end

fid = fopen('C:\Work\ICR_LCM\BR1458_smallRemoved.rois', 'w');
roiList.outputXML(fid, 0);
fclose(fid);



%%

hande = imread('C:\Work\ICR_LCM\ICR\BR1458_5_H&E.jpg');

cpselect(hande, k5./max(k5));

%%

optical = imread('C:\Work\ICR_LCM\ICR\BR1458_5_zoom.JPG');

cpselect(optical, k5./max(k5));

%%

k5rgb = ind2rgb(k5+1, parula(6));

figure, imagesc(k5rgb);
axis image; axis off;
%%

transform = fitgeotrans(fixedPoints, movingPoints, 'affine');

k5_optical = imwarp(k5rgb, transform, 'OutputView', imref2d(size(optical)));

fusedImage = imfuse(imresize(k5_optical, 0.25), imresize(optical, 0.25), 'blend');

figure, %imagesc(k5_optical);
imagesc(fusedImage);
axis image; axis off;

%%

registration = [];

registration.fixedImage = 'C:\Work\ICR_LCM\ICR\BR1458_5_zoom.JPG';

for i = 1:size(movingPoints, 1)
    registration.fixedPoints(i).x = movingPoints(i, 1);
    registration.fixedPoints(i).y = movingPoints(i, 2);
end
registration.movingImage = 'C:\Work\ICR_LCM\ICR\data\BR1458_5_cluster\';

for i = 1:size(fixedPoints, 1)
    registration.movingPoints(i).x = fixedPoints(i, 1);
    registration.movingPoints(i).y = fixedPoints(i, 2);
end

fid = fopen('C:\Work\ICR_LCM\ICR\BR1458_registration.regi', 'w');
fwrite(fid, jsonencode(registration));
fclose(fid);
%%
a = [fixedPoints ones(size(fixedPoints, 1), 1)]';
b = [movingPoints ones(size(movingPoints, 1), 1)]';

movingToFixed =  b' \ a'

movingToFixed' * [movingPoints(1, :) 1]'