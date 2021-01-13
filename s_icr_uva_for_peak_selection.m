
path = 'X:\ICR Breast PDX\Data Processing Outputs\pdx & primary tumours\uva\neg desi 3d pdx & primary\tissue only\no norm\';

%%

load([path ' wt vs mutant (3d pdx only 4 pairs)'])

data = ones(size(table,1)-1,4);

disp(table(1,1));
data(:,1) = double(table(2:end,1)); % AUC (t-1282-4 vs t-1458-2)

disp(table(1,6));
data(:,2) = double(table(2:end,6)); % AUC (t-1282-5 vs t-1458-4)

disp(table(1,11));
data(:,3) = double(table(2:end,11)); % AUC (t-1282-3 vs t-1458-3)

disp(table(1,16));
data(:,4) = double(table(2:end,16)); % AUC (t-1282-2 vs t-1458-5)

disp(table(1,21));
mzs = double(table(2:end,21));

logical_data_wt = data>=0.7;
logical_data_mt = data<=0.3;

logical_data = logical_data_wt + logical_data_mt;

mzs4mva = unique(mzs(sum(logical_data,2)==4,1)); % 306 peaks

save([path 'mzs4mva_wtvsmt_bothways_across4pairs.mat'],'mzs4mva')

% mzs4mva_3dpdx = unique(mzs(sum(logical_data,2)==4,1)); % 311 peaks

load('X:\ICR Breast PDX\Data Processing Outputs\pdx & primary tumours\spectra details\2019_10_09_PDX_BR1458_2_NEG_100ss Analyte 3\tissue only\datacubeonly_peakDetails')

A = repmat(mzs4mva', size(datacubeonly_peakDetails,1), 1);
B = datacubeonly_peakDetails(:,2);
[ minValue, closestIndex ] = min(abs(A-B));

ppm = minValue./mzs4mva'.*1e6;

mzs4mva_3dpdxandprimary = datacubeonly_peakDetails(closestIndex,2);
mzs4mva_3dpdxandprimary = mzs4mva_3dpdxandprimary(ppm<=30);
 
% save([path 'mzs4mva_wtvsmt_bothways_across4pairs.mat'],'mzs4mva_3dpdxandprimary') % 278 peaks

%%

load([path 'cluster vs cluster (3d pdx only kmeans 8)'])

data = ones(size(table,1)-1,2);

disp(table(1,51));
data(:,1) = double(table(2:end,51)); % AUC (c5 vs c2)

disp(table(1,146));
data(:,2) = double(table(2:end,146)); % AUC (c2 vs c5)

disp(table(1,281));
mzs = double(table(2:end,281));

logical_data_wt = data>=0.7;
logical_data_mt = data<=0.3;

logical_data = logical_data_wt + logical_data_mt;

mzs4mva = unique(mzs(sum(logical_data,2)>=1,1));

save([path 'mzs4mva_c5vsc2_bothways.mat'],'mzs4mva') % 447 peaks
