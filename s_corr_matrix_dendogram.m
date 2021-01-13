
% Load Clustering outputs

% load('X:\ICR Breast PDX\Data Processing Outputs\3d pdx pilot\mva 4000 highest peaks\neg desi 3d pdx pilot\tissue only\kmeans 8 components\no norm\C')
% load('X:\ICR Breast PDX\Data Processing Outputs\3d pdx pilot\mva 4000 highest peaks\neg desi 3d pdx pilot\tissue only\kmeans 8 components\no norm\datacube_mzvalues_indexes')
% load('X:\ICR Breast PDX\Data Processing Outputs\3d pdx pilot\spectra details\2019_10_09_PDX_BR1458_2_NEG_100ss Analyte 3\tissue only\datacubeonly_peakDetails')

load('X:\ICR Breast PDX\Data Processing Outputs\pdx & primary tumours\mva 4000 highest peaks\neg desi 3d pdx & primary\tissue only\kmeans 20 components\no norm\C')
load('X:\ICR Breast PDX\Data Processing Outputs\pdx & primary tumours\mva 4000 highest peaks\neg desi 3d pdx & primary\tissue only\kmeans 20 components\no norm\datacube_mzvalues_indexes')
load('X:\ICR Breast PDX\Data Processing Outputs\pdx & primary tumours\spectra details\2019_10_09_PDX_BR1458_2_NEG_100ss Analyte 3\tissue only\datacubeonly_peakDetails')

load('X:\ICR Breast PDX\Data Processing Outputs\pdx & primary tumours\mva 4000 highest peaks\neg desi 3d pdx & primary\tissue only\kmeans 20 components\tic\C')
load('X:\ICR Breast PDX\Data Processing Outputs\pdx & primary tumours\mva 4000 highest peaks\neg desi 3d pdx & primary\tissue only\kmeans 20 components\tic\datacube_mzvalues_indexes')
load('X:\ICR Breast PDX\Data Processing Outputs\pdx & primary tumours\spectra details\2019_10_09_PDX_BR1458_2_NEG_100ss Analyte 3\tissue only\datacubeonly_peakDetails')

load('X:\ICR Breast PDX\Data Processing Outputs\pdx & primary tumours\mva 4000 highest peaks\neg desi 3d pdx & primary\tissue only\kmeans 20 components\RMS\C')
load('X:\ICR Breast PDX\Data Processing Outputs\pdx & primary tumours\mva 4000 highest peaks\neg desi 3d pdx & primary\tissue only\kmeans 20 components\RMS\datacube_mzvalues_indexes')
load('X:\ICR Breast PDX\Data Processing Outputs\pdx & primary tumours\spectra details\2019_10_09_PDX_BR1458_2_NEG_100ss Analyte 3\tissue only\datacubeonly_peakDetails')


mz = datacubeonly_peakDetails(datacube_mzvalues_indexes,2);
spectra = C';

cmap = tab20(20);
figure; for k = 1:20; hold on; stem(mz,spectra(:,k),'d','filled','color',cmap(k,:)); legend; end

%%

% Correlation matrix

corr_matrix = corr(spectra);
cmap = gray;
figure;
h = heatmap(round(corr_matrix,2),'Colormap',cmap); caxis([0 1]); 

h.Title = 'Spectral Correlation';
h.XLabel = 'Cluster ID';
h.YLabel = 'Cluster ID';

%%

% Dendrogram plots

% Warning: Non-monotonic cluster tree -- the median linkage is probably not appropriate.
% Warning: Non-monotonic cluster tree -- the centroid linkage is probably not appropriate.

figure; 
% subplot(3,4,1); h = dendrogram(linkage(C,'single','euclidean')); title('single - euclidean'); grid on;
% subplot(3,4,2); h = dendrogram(linkage(C,'complete','euclidean')); title('complete - euclidean'); grid on;
subplot(3,4,1); h = dendrogram(linkage(C,'average','euclidean')); title('average - euclidean'); grid on;
subplot(3,4,2); h = dendrogram(linkage(C,'weighted','euclidean')); title('weighted - euclidean'); grid on;
% subplot(3,5,4); h = dendrogram(linkage(C,'median','euclidean')); title('median - euclidean'); grid on;
subplot(3,4,3); h = dendrogram(linkage(C,'centroid','euclidean')); title('centroid - euclidean'); grid on;
subplot(3,4,4); h = dendrogram(linkage(C,'ward','euclidean')); title('ward - euclidean'); grid on;

% subplot(3,4,1+4); h = dendrogram(linkage(C,'single','cosine')); title('single - cosine'); grid on;
% subplot(3,4,2+4); h = dendrogram(linkage(C,'complete','cosine')); title('complete - cosine'); grid on;
subplot(3,4,1+4); h = dendrogram(linkage(C,'average','cosine')); title('average - cosine'); grid on;
subplot(3,4,2+4); h = dendrogram(linkage(C,'weighted','cosine')); title('weighted - cosine'); grid on;
% subplot(3,5,4+5); h = dendrogram(linkage(C,'median','cosine')); title('median - cosine'); grid on;
% subplot(3,5,5+5); h = dendrogram(linkage(C,'ward','cosine')); title('ward - cosine'); grid on;

% subplot(3,4,1+8); h = dendrogram(linkage(C,'single','correlation')); title('single - correlation'); grid on;
% subplot(3,4,2+8); h = dendrogram(linkage(C,'complete','correlation')); title('complete - correlation'); grid on;
subplot(3,4,1+8); h = dendrogram(linkage(C,'average','correlation')); title('average - correlation'); grid on;
subplot(3,4,2+8); h = dendrogram(linkage(C,'weighted','correlation')); title('weighted - correlation'); grid on;
% subplot(3,5,4+10); h = dendrogram(linkage(C,'median','correlation')); title('median - correlation'); grid on;
% subplot(3,5,5+10); h = dendrogram(linkage(C,'ward','correlation')); title('ward - correlation'); grid on;

figure;
hidx = cluster(linkage(C,'weighted','euclidean'),'criterion','distance','maxclust',4);

figure;
subplot(1,4,1)
linkage_output = linkage(C,'single','cosine');
[ H, T, outperm ] = dendrogram(linkage_output, 'Orientation', 'left');
set(H,'color',[0 0 0])
title('Dendogram (single linkage, cosine distance)'); 
grid on;
subplot(1,4,2:4)
file = 'X:\ICR Breast PDX\Data Processing Outputs\pdx & primary tumours\mva 4000 highest peaks\neg desi 3d pdx & primary\tissue only\kmeans 20 components\no norm\clusters_table.txt';
txt_row = strcat(repmat('%s\t',1,size(table,2)-1),'%s\n');        
fileID = fopen(file,'r');
table = fscanf(fileID,'%s\t');
fclose(fileID);
for row = 1:size(clusters_vs_samples,1)

    clusterids = clusters_vs_samples.clusterId
    
end


%%

K = 8;
load('X:\ICR Breast PDX\Data Processing Outputs\3d pdx pilot\uva\neg desi 3d pdx\tissue only\no norm\cluster vs cluster (kmeans 8)')

uva_output = table;

auc_data0 = double(uva_output(2:end, [1:5:5*(K*K-K) find(strcmp(uva_output(1,:),'meas mz'))]));
[~, rows] = unique(uva_output(:,strcmp(uva_output(1,:),'meas mz')));

rows = rows(2:end);

M07 = NaN*ones(K,K);
M03 = NaN*ones(K,K);
for col = 1:5:5*(K*K-K)
    B0 = regexp(uva_output(1, col),'\d*','Match');
    B = double(B0);
    M07(B(1),B(2)) = sum(double(uva_output(rows,col))>=0.7);
    M03(B(1),B(2)) = sum(double(uva_output(rows,col))<=0.3);
end

cmap = gray;
figure;
subplot(1,2,1); h = heatmap(M07,'Colormap',cmap(end:-1:1,:)); caxis([0 size(rows,1)]); 
h.Title = '# peaks with AUC>=0.7 (Row Index vs Col Index)';
h.XLabel = 'Cluster ID';
h.YLabel = 'Cluster ID';
subplot(1,2,2); h = heatmap(M07+M03,'Colormap',cmap(end:-1:1,:)); caxis([0 size(rows,1)]); 
h.Title = '# peaks with AUC>=0.7 or AUC<=0.3';
h.XLabel = 'Cluster ID';
h.YLabel = 'Cluster ID';

%% Other distances using the groups of pixels within the clusters instead of the mean cluster

K = 8;
load('X:\ICR Breast PDX\Data Processing Outputs\3d pdx pilot\similarity analysis\neg desi 3d pdx pilot\tissue only\no norm\cluster vs cluster (kmeans 8)')

sa_output = table;

names1 = sa_output(1, 1:2:2*(K*K-K));
names2 = sa_output(1, 2:2:2*(K*K-K));

dist1 = double(sa_output(2, 1:2:2*(K*K-K)));
dist2 = double(sa_output(2, 2:2:2*(K*K-K)));

M1 = zeros(K,K);
M2 = zeros(K,K);
for col = 1:size(dist1,2)
    B0 = regexp(names1(1, col),'\d*','Match');
    B = double(B0);
    M1(B(1),B(2)) = dist1(1,col);
end
for col = 1:size(dist2,2)
    B0 = regexp(names2(1, col),'\d*','Match');
    B = double(B0);
    M2(B(1),B(2)) = dist2(1,col);
end

M1 = round(mat2gray(M1),2);
M2 = round(mat2gray(M2),2);

cmap = gray;
figure;
subplot(1,2,1); h = heatmap(M1,'Colormap',cmap(end:-1:1,:)); %caxis([0 size(rows,1)]); 
h.Title = 'Wasserstein Dist';
h.XLabel = 'Cluster ID';
h.YLabel = 'Cluster ID';
subplot(1,2,2); h = heatmap(M2,'Colormap',cmap(end:-1:1,:)); %caxis([0 size(rows,1)]); 
h.Title = 'Xaviers Wasserstein Dist';
h.XLabel = 'Cluster ID';
h.YLabel = 'Cluster ID';


    
