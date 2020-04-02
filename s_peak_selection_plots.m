
datacubeonly_peakDetails99 = datacubeonly_peakDetails(datacube_indexes_99,:);
datacubeonly_peakDetailsnone = datacubeonly_peakDetails(datacube_indexes_none,:);

%%

samples_num = 0;
max_position_vector = [];

for i = 1:100
    y = totalSpectrum_intensities(datacubeonly_peakDetails99(i,6):datacubeonly_peakDetails99(i,7));
    [ ~, max_position ] = max(y);
    max_position_vector(i,1) = max_position;
    samples_num = max(samples_num, size(y,2));
end

new_y = NaN*ones(100,2*samples_num);
for i = 1:100
    y = totalSpectrum_intensities(datacubeonly_peakDetails99(i,6):datacubeonly_peakDetails99(i,7));    
    new_y(i,[1:size(y,2)] + samples_num - max_position_vector(i,1)) = y;
end

norm_y = (new_y-min(new_y,[],2))./max(new_y-min(new_y,[],2),[],2);

figure; plot(1:356,norm_y)

%%

samples_num = 0;
max_position_vector = [];

for i = 1:100
    y = totalSpectrum_intensities(datacubeonly_peakDetailsnone(i,6):datacubeonly_peakDetailsnone(i,7));
    [ ~, max_position ] = max(y);
    max_position_vector(i,1) = max_position;
    samples_num = max(samples_num, size(y,2));
end

new_y = NaN*ones(100,2*samples_num);
for i = 1:100
    y = totalSpectrum_intensities(datacubeonly_peakDetailsnone(i,6):datacubeonly_peakDetailsnone(i,7));    
    new_y(i,[1:size(y,2)] + samples_num - max_position_vector(i,1)) = y;
end

norm_y = (new_y-min(new_y,[],2))./max(new_y-min(new_y,[],2),[],2);

figure; plot(1:356,norm_y)

%%

load('X:\Beatson\dpo mva article\negative DESI\mva 100 highest peaks\negative DESI small intestine\tissue only\tsne\pqn median\idx.mat')

idx_100 = idx;

cluster_100_apckras = idx;
cluster_100_apckras(idx==0) = [];
cluster_100_apckras(cluster_100_apckras~=3) = 0;
cluster_100_apckras(cluster_100_apckras==3) = 1;

cluster_100_apc = idx;
cluster_100_apc(idx==0) = [];
cluster_100_apc(cluster_100_apc~=8) = 0;
cluster_100_apc(cluster_100_apc==8) = 1;

%%

load('X:\Beatson\dpo mva article\negative DESI\mva Shorter Beatson metabolomics & CRUK list\negative DESI small intestine\tissue only\tsne\pqn median\idx.mat')

idx_sl = idx;

cluster_sl_apckras = idx;
cluster_sl_apckras(idx==0) = [];
cluster_sl_apckras(cluster_sl_apckras~=14) = 0;
cluster_sl_apckras(cluster_sl_apckras==14) = 1;

cluster_sl_apc = idx;
cluster_sl_apc(idx==0) = [];
cluster_sl_apc(cluster_sl_apc~=1) = 0;
cluster_sl_apc(cluster_sl_apc==1) = 1;

%%

similarity = dice(cluster_100_apckras,cluster_sl_apckras)
similarity = dice(cluster_100_apc,cluster_sl_apc)


dice(idx_100,idx_sl)
