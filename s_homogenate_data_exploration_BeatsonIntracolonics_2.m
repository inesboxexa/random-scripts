
% This script was written by Teresa Murta for the
% 2020110334-NiCEACR measurement service job.

% % This script was written to explore the homogenate data collected before
% % and after most of the sample of interest data, with the aim to understand
% % the relationship between the drug intensity trends seen in the samples of
% % interest with the intensity of severel ions in the homogenate.
%
% % It runs the following short data analysis pipeline:
% % (1) definition of the groups of samples that need to be combined (e.g.:
% % all time points for each patient)
% % (2) for each group of samples:
% % (2.1) generating the SIMS-ToF total spectrum of each dataset within the
% % group, and summing them up to generate the representative total spectrum
% % for that group of samples (this step guaranties the same mass axis for
% % all SIMS-ToF datasets that need to be explored)
% % (2.2) peak picking the representative SIMS-ToF total spectrum of a given
% % group of samples
% % (2.3) generating the ion images for the most intense 500 ions, and
% % gathering them in a unique matlab variable
% % (3) look at the mean, median, percentile 25 and 75 of a list of 12 ions
% % (defined by the experimentalist) across the samples (e.g. time points)
% % (4) running a principal componenet analysis (PCA) using the intensities
% % of the 500 most intense ions to visualise the trends across the samples

%%

normalisation = "tic";

%% Gathering the veal-brain-homogenate and tumours data

data_type = 'vbh & tumours';

data_folders = { 'D:\veal-brain-homogenate-study\vbh data\' }; dataset_name_portion = { '*' };
filesToProcess = []; for i = 1:length(data_folders); filesToProcess = [ filesToProcess; dir([data_folders{i} filesep dataset_name_portion{i} '.imzML']) ]; end % Files and adducts information gathering
MyFilesInfo1 = filesToProcess([ 16 21 33 ]);

data_folders = { 'D:\veal-brain-homogenate-study\tumours data\'}; dataset_name_portion = { '*' };
filesToProcess = []; for i = 1:length(data_folders); filesToProcess = [ filesToProcess; dir([data_folders{i} filesep dataset_name_portion{i} '.imzML']) ]; end % Files and adducts information gathering
MyFilesInfo2 = filesToProcess([ 10:11 14:20 ]);

MyFilesInfo = [
    MyFilesInfo1(1);
    MyFilesInfo2(1:2);
    MyFilesInfo1(2);
    MyFilesInfo2(3:5);
    MyFilesInfo1(3);
    MyFilesInfo2(6:end);
    ];

gemm_id = [];
date_id = [];
replicate_id = [];
dataMatrix = [];
sample_ids = [];
for i = 1:size(MyFilesInfo,1)
    
    file_name = strrep(MyFilesInfo(i).name(1:end-6), '_','-');
    
    % Spectral details
    
    load([ 'D:\veal-brain-homogenate-study\dpo\spectra details' filesep MyFilesInfo(i).name(1:end-6) '\tissue only\' char(normalisation) filesep 'data' ]);
    load([ 'D:\veal-brain-homogenate-study\dpo\spectra details' filesep MyFilesInfo(i).name(1:end-6) '\tissue only\datacubeonly_peakDetails' ]);
    load([ 'D:\veal-brain-homogenate-study\dpo\spectra details' filesep MyFilesInfo(i).name(1:end-6) '\tissue only\width' ]);
    load([ 'D:\veal-brain-homogenate-study\dpo\spectra details' filesep MyFilesInfo(i).name(1:end-6) '\tissue only\height' ]);
    load([ 'D:\veal-brain-homogenate-study\dpo\spectra details' filesep MyFilesInfo(i).name(1:end-6) '\tissue only\totalSpectrum_intensities' ]);
    load([ 'D:\veal-brain-homogenate-study\dpo\spectra details' filesep MyFilesInfo(i).name(1:end-6) '\tissue only\totalSpectrum_mzvalues' ]);
    
    % ROIs
    
    switch MyFilesInfo(i).name(1:end-6)
        
        case '20201127_DESI_NEG_VBH_2x2_3'
            load('D:\veal-brain-homogenate-study\dpo\rois\20201127_DESI_NEG_VBH_2x2_3\vbh-20201127\roi')
            mask = reshape(roi.pixelSelection',[],1);
            roi_masks = mask;
            gemm_id = [gemm_id, 1];
            date_id = [date_id, 1];
            replicate_id = [replicate_id, 1];
            
        case '20201127_intracolonics_neg_xDESI_R2B2_s19_1'
            load('D:\veal-brain-homogenate-study\dpo\rois\20201127_intracolonics_neg_xDESI_R2B2_s19_1\apc-1-R2-B2-S19\roi')
            mask = reshape(roi.pixelSelection',[],1);
            roi_masks = mask;
            gemm_id = [gemm_id, 2];
            date_id = [date_id, 1];
            replicate_id = [replicate_id, 2];
            load('D:\veal-brain-homogenate-study\dpo\rois\20201127_intracolonics_neg_xDESI_R2B2_s19_1\apc-kras-2-R2-B2-S19\roi')
            mask = reshape(roi.pixelSelection',[],1);
            roi_masks = [roi_masks, mask];
            gemm_id = [gemm_id, 3];
            date_id = [date_id, 1];
            replicate_id = [replicate_id, 4];
            
        case '20201127_intracolonics_neg_xDESI_R2B2_s19_2'
            load('D:\veal-brain-homogenate-study\dpo\rois\20201127_intracolonics_neg_xDESI_R2B2_s19_2\apc-4-R2-B2-S19\roi')
            mask = reshape(roi.pixelSelection',[],1);
            roi_masks = mask;
            gemm_id = [gemm_id, 2];
            date_id = [date_id, 1];
            replicate_id = [replicate_id, 3];
            load('D:\veal-brain-homogenate-study\dpo\rois\20201127_intracolonics_neg_xDESI_R2B2_s19_2\apc-kras-3-R2-B2-S19\roi')
            mask = reshape(roi.pixelSelection',[],1);
            roi_masks = [roi_masks, mask];
            gemm_id = [gemm_id, 3];
            date_id = [date_id, 1];
            replicate_id = [replicate_id, 5];
            load('D:\veal-brain-homogenate-study\dpo\rois\20201127_intracolonics_neg_xDESI_R2B2_s19_2\apc-kras-5-R2-B2-S19\roi')
            mask = reshape(roi.pixelSelection',[],1);
            roi_masks = [roi_masks, mask];
            gemm_id = [gemm_id, 3];
            date_id = [date_id, 1];
            replicate_id = [replicate_id, 6];
            load('D:\veal-brain-homogenate-study\dpo\rois\20201127_intracolonics_neg_xDESI_R2B2_s19_2\apc-kras-6-R2-B2-S19\roi')
            mask = reshape(roi.pixelSelection',[],1);
            roi_masks = [roi_masks, mask];
            gemm_id = [gemm_id, 3];
            date_id = [date_id, 1];
            replicate_id = [replicate_id, 7];
            
        case '20201202_DESI_NEG_VBH_2x2_3'
            load('D:\veal-brain-homogenate-study\dpo\rois\20201202_DESI_NEG_VBH_2x2_3\vbh-20201202\roi')
            mask = reshape(roi.pixelSelection',[],1);
            roi_masks = mask;
            gemm_id = [gemm_id, 1];
            date_id = [date_id, 2];
            replicate_id = [replicate_id, 1];
            
        case '20201202_intracolonics_neg_xDESI_R2B2_s20_1'
            load('D:\veal-brain-homogenate-study\dpo\rois\20201202_intracolonics_neg_xDESI_R2B2_s20_1\apc-1-R2-B2-S20\roi')
            mask = reshape(roi.pixelSelection',[],1);
            roi_masks = mask;
            gemm_id = [gemm_id, 2];
            date_id = [date_id, 2];
            replicate_id = [replicate_id, 2];
            load('D:\veal-brain-homogenate-study\dpo\rois\20201202_intracolonics_neg_xDESI_R2B2_s20_1\apc-kras-2-R2-B2-S20\roi')
            mask = reshape(roi.pixelSelection',[],1);
            roi_masks = [roi_masks, mask];
            gemm_id = [gemm_id, 3];
            date_id = [date_id, 2];
            replicate_id = [replicate_id, 4];
            
        case '20201202_intracolonics_neg_xDESI_R2B2_s20_2'
            load('D:\veal-brain-homogenate-study\dpo\rois\20201202_intracolonics_neg_xDESI_R2B2_s20_2\apc-4-R2-B2-S20\roi')
            mask = reshape(roi.pixelSelection',[],1);
            roi_masks = mask;
            gemm_id = [gemm_id, 2];
            date_id = [date_id, 2];
            replicate_id = [replicate_id, 3];
            load('D:\veal-brain-homogenate-study\dpo\rois\20201202_intracolonics_neg_xDESI_R2B2_s20_2\apc-kras-3-R2-B2-S20\roi')
            mask = reshape(roi.pixelSelection',[],1);
            roi_masks = [roi_masks, mask];
            gemm_id = [gemm_id, 3];
            date_id = [date_id, 2];
            replicate_id = [replicate_id, 5];
            
        case '20201202_intracolonics_neg_xDESI_R2B2_s20_3'
            load('D:\veal-brain-homogenate-study\dpo\rois\20201202_intracolonics_neg_xDESI_R2B2_s20_3\apc-kras-5-R2-B2-S20\roi')
            mask = reshape(roi.pixelSelection',[],1);
            roi_masks = mask;
            gemm_id = [gemm_id, 3];
            date_id = [date_id, 2];
            replicate_id = [replicate_id, 6];
            
        case '20201207_xDESI_V3_NEG_VBH'
            load('D:\veal-brain-homogenate-study\dpo\rois\20201207_xDESI_V3_NEG_VBH\vbh-20201207\roi')
            mask = reshape(roi.pixelSelection',[],1);
            roi_masks = mask;
            gemm_id = [gemm_id, 1];
            date_id = [date_id, 3];
            replicate_id = [replicate_id, 1];
            
        case '20201207_intracolonics_neg_xDESI_R2B2_s21_t1_2'
            load('D:\veal-brain-homogenate-study\dpo\rois\20201207_intracolonics_neg_xDESI_R2B2_s21_t1_2\apc-1-R2-B2-S21\roi')
            mask = reshape(roi.pixelSelection',[],1);
            roi_masks = mask;
            gemm_id = [gemm_id, 2];
            date_id = [date_id, 3];
            replicate_id = [replicate_id, 2];
            load('D:\veal-brain-homogenate-study\dpo\rois\20201207_intracolonics_neg_xDESI_R2B2_s21_t1_2\apc-kras-2-R2-B2-S21\roi')
            mask = reshape(roi.pixelSelection',[],1);
            roi_masks = [roi_masks, mask];
            gemm_id = [gemm_id, 3];
            date_id = [date_id, 3];
            replicate_id = [replicate_id, 4];
            
        case '20201207_intracolonics_neg_xDESI_R2B2_s21_t3_4'
            load('D:\veal-brain-homogenate-study\dpo\rois\20201207_intracolonics_neg_xDESI_R2B2_s21_t3_4\apc-4-R2-B2-S21\roi')
            mask = reshape(roi.pixelSelection',[],1);
            roi_masks = mask;
            gemm_id = [gemm_id, 2];
            date_id = [date_id, 3];
            replicate_id = [replicate_id, 3];
            load('D:\veal-brain-homogenate-study\dpo\rois\20201207_intracolonics_neg_xDESI_R2B2_s21_t3_4\apc-kras-3-R2-B2-S21\roi')
            mask = reshape(roi.pixelSelection',[],1);
            roi_masks = [roi_masks, mask];
            gemm_id = [gemm_id, 3];
            date_id = [date_id, 3];
            replicate_id = [replicate_id, 5];
            
        case '20201207_intracolonics_neg_xDESI_R2B2_s21_t5'
            load('D:\veal-brain-homogenate-study\dpo\rois\20201207_intracolonics_neg_xDESI_R2B2_s21_t5\apc-kras-5-R2-B2-S21\roi')
            mask = reshape(roi.pixelSelection',[],1);
            roi_masks = mask;
            gemm_id = [gemm_id, 3];
            date_id = [date_id, 3];
            replicate_id = [replicate_id, 6];
            
        case '20201207_intracolonics_neg_xDESI_R2B2_s21_t6'
            load('D:\veal-brain-homogenate-study\dpo\rois\20201207_intracolonics_neg_xDESI_R2B2_s21_t6\apc-kras-6-R2-B2-S21\roi')
            mask = reshape(roi.pixelSelection',[],1);
            roi_masks = mask;
            gemm_id = [gemm_id, 3];
            date_id = [date_id, 3];
            replicate_id = [replicate_id, 7];
            
    end
    
    for roii = 1:size(roi_masks,2)
        
        dataMatrix0 = data.*roi_masks(:,roii);
        dataMatrix0(dataMatrix0==0) = NaN;
        mzs0 = datacubeonly_peakDetails(:,2);
        
        if isempty(dataMatrix)
            dataMatrix = dataMatrix0;
            mzs = mzs0;
            spectraMatrix = totalSpectrum_intensities';
            spectra_mzs = totalSpectrum_mzvalues';
        else
            if size(dataMatrix,1)<size(dataMatrix0,1)
                dataMatrix = cat(1,dataMatrix,NaN*ones(size(dataMatrix0,1)-size(dataMatrix,1),size(dataMatrix,2),size(dataMatrix,3)));
            elseif size(dataMatrix,1)>size(dataMatrix0,1)
                dataMatrix0 = cat(1,dataMatrix0,NaN*ones(size(dataMatrix,1)-size(dataMatrix0,1),size(dataMatrix0,2),size(dataMatrix0,3)));
            end
            dataMatrix = cat(3,dataMatrix,dataMatrix0);
            mzs = [mzs,mzs0];
            spectraMatrix = [spectraMatrix,totalSpectrum_intensities'];
        end
        
        sample_ids = [ sample_ids, string(file_name) ];
        
    end
    
end

%% Ploting the mean and median intensities of the 12 peaks selected as reference peaks by the experimentalist

peaksOfInterest_cruk = [
    174.04
    89.024
    146.046
    135.03
    132.03
    145.062
    218.102
    133.014
    111.02
    179.056
    267.072
    124.006
    ];

peaksOfInterest_lipids = [
    303.232
    327.232
    281.248
    329.248
    331.264
    788.542
    788.542
    305.246
    283.264
    255.232
    834.526
    885.546
    ];

peaksOfInterest_pc2 = [
    306.076
    482.974
    630.492
    645.448
    673.536
    679.504
    734.56
    ];

peaksOfInterest_lowestp_CRUK = [
    500.278
    202.11
    306.076
    317.174
    303.05
    145.014
    385.348
    322.18
    448.302
    427.25
    303.232
    239.058
    411.278
    210.016
    134.046
    267.072
    120.012
    130.086
    164.072
    279.232
    407.28
    886.08
    307
    173.008
    180.066
    252.094
    124.006
    484.272
    426.034
    373.042
    266.114
    117.018
    282.084
    303.082
    263.016
    461.362
    242.078
    318.06
    347.046
    236.018
    174.04
    315.208
    278.054
    433.1
    227.036
    203.082
    339.208
    346.054
    226.996
    135.03
    438.966
    208.982
    440.134
    289.036
    289.116
    325.008
    114.056
    103.002
    353.162
    337.06
    139.076
    383.008
    255.232
    181.07
    401.022
    209.082
    179.056
    254.082
    808.132
    210.06
    663.1
    167.02
    200.054
    946.144
    186.018
    516.96
    436.996
    342.976
    154.062
    367.008
    191.02
    183.046
    184.986
    339.07
    166.976
    166.976
    132.03
    138.978
    243.062
    257.03
    146.046
    381.134
    291.204
    279.038
    302.066
    521.988
    426.01
    844.108
    134.982
    239.016
    180.99
    111.02
    159.984
    145.062
    272.004
    221.06
    489.03
    174.088
    808.112
    946.128
    461.99
    202.946
    ];
peaksOfInterest_lowestp = [
    500.278
    500.278
    500.278
    500.278
    500.278
    235.168
    235.168
    244.136
    244.136
    244.136
    244.136
    290.1
    290.1
    679.504
    679.504
    679.504
    582.28
    329.192
    677.49
    677.49
    677.49
    690.532
    690.532
    528.308
    528.308
    356.262
    897.752
    897.752
    717.562
    717.562
    717.562
    717.562
    731.114
    731.114
    973.782
    973.782
    629.49
    629.49
    629.49
    629.49
    629.49
    949.782
    949.782
    202.11
    202.11
    202.11
    273.184
    733.558
    733.558
    733.558
    733.558
    733.558
    470.318
    470.318
    469.314
    469.314
    899.766
    501.28
    501.28
    678.492
    678.492
    543.298
    733.112
    478.254
    478.254
    443.278
    443.278
    921.752
    ];

peaksOfInterest_specialList = [
    145.062
    181.036
    132.03
    146.046
    179.056
    215.032
    ];

%%

peaksOfInteresti = 0;
for peak = peaksOfInterest_lipids'
    [~, index] = ismembertol(mzs0,peak,1e-8);
    disp(sum(index))
    peaksOfInteresti = peaksOfInteresti+index;
end

peaksOfInteresti = unique(find(peaksOfInteresti));
peaksOfInteresti = peaksOfInteresti(1:6);

% Update data and spectral channels

dataMatrix2plot = dataMatrix(:,peaksOfInteresti,:);
% dataMatrix2plot = dataMatrix(:,152,:)./dataMatrix(:,112,:); % Ratio
mzs2plot = mzs0(peaksOfInteresti);
% mzs2plot = mzs0(152)/mzs0(112); % Ratio

means = squeeze(mean(dataMatrix2plot(:,:,:),1,'omitnan'));
medians = squeeze(median(dataMatrix2plot(:,:,:),1,'omitnan'));
prctile75 = squeeze(prctile(dataMatrix2plot(:,:,:),75,1));
prctile25 = squeeze(prctile(dataMatrix2plot(:,:,:),25,1));

% means = squeeze(mean(dataMatrix2plot(:,:,:),1,'omitnan'))'; % Ratio
% medians = squeeze(median(dataMatrix2plot(:,:,:),1,'omitnan'))'; % Ratio
% prctile75 = squeeze(prctile(dataMatrix2plot(:,:,:),75,1))'; % Ratio
% prctile25 = squeeze(prctile(dataMatrix2plot(:,:,:),25,1))'; % Ratio

colors = viridis(5);
colors(1,:) = [0 0 0];

shapes={'d','o','s'};
figure
for peaki = 1:size(means,1)
    subplot(3,2,peaki)
    hold on
    title(['m/z ',  num2str(round(mzs2plot(peaki),3))])
    xi=0;
    for replicatei = unique(replicate_id)
        for gemmi = unique(gemm_id)
            for datei = unique(date_id)
                columni = logical((replicatei==replicate_id) .* (gemmi==gemm_id) .* (datei==date_id));
                if sum(columni)==1
                    xi=xi+1;
                    plot(xi,means(peaki,columni),'x','color', colors(gemmi,:),'MarkerEdgeColor',colors(gemmi,:),'MarkerFaceColor',[1 1 1])
                    plot(xi,medians(peaki,columni),shapes{datei},'color', colors(gemmi,:),'MarkerEdgeColor',colors(gemmi,:),'MarkerFaceColor',colors(gemmi,:))
                    plot(xi,prctile25(peaki,columni),'^','color', colors(gemmi,:),'MarkerEdgeColor',colors(gemmi,:),'MarkerFaceColor',[1 1 1])
                    plot(xi,prctile75(peaki,columni),'v','color', colors(gemmi,:),'MarkerEdgeColor',colors(gemmi,:),'MarkerFaceColor',[1 1 1])
                end
            end
        end
    end
    
    stem([3.5:3:20], 2*max(means(peaki,:))*ones(1,6),'k')
    axis([0 length(replicate_id)+1 0 1.1*max(prctile75(peaki,:))])
    
    grid on
    ylabel('ion counts')
    xlabel('sample')
    xticks(1:xi)
    switch data_type
        case 'vbh & tumours'
            xticklabels({...
                '27-Nov','2-Dec','7-Dec',...
                '27-Nov','2-Dec','7-Dec',...
                '27-Nov','2-Dec','7-Dec',...
                '27-Nov','2-Dec','7-Dec',...
                '27-Nov','2-Dec','7-Dec',...
                '27-Nov','2-Dec','7-Dec',...
                '27-Nov','7-Dec',...
                });
    end
    xtickangle(45)
    if peaki == 2 % size(totalcount,1)
        legend({'mean','median','percentile 25','percentile 75'})
        legend('boxoff')
        legend('Location','best')
    end
end

%% Running a Principal Component Analysis (PCA) of the 500 most intense ions and plotting the scores which illustrate the intensaty trends across samples

colors = tab10;
figure
data4pca = squeeze(mean(dataMatrix,1,'omitnan'))';
data4pca(:,sum(isnan(data4pca))>0) = [];
data4pca = (data4pca-mean(data4pca,'omitnan'))./std(data4pca(sum(~isnan(data4pca)==0,2)==0,:));
[coeff,score,latent,tsquared,explained,mu] = pca(data4pca);
bar(score(:,1:5),'EdgeColor',[1 1 1])
hold on
stem([3.5:3:20], 2*max(abs(score(:)))*ones(1,6),'k')
stem([3.5:3:20], -2*max(abs(score(:)))*ones(1,6),'k')
axis([0 length(replicate_id)+1 1.1*min(score(:)) 1.1*max(score(:))])
legend({['PC 1 (', num2str(round(explained(1),2)) '% var expl)'], ['PC 2 (', num2str(round(explained(2),2)) '% var expl)'], ['PC 3 (', num2str(round(explained(3),2)) '% var expl)'], ['PC 4 (', num2str(round(explained(4),2)) '% var expl)'], ['PC 5 (', num2str(round(explained(5),2)) '% var expl)']});
legend('boxoff')
legend('Location','best')
xticks(1:size(dataMatrix,3))
switch data_type
    case 'vbh & tumours'
        xticklabels({...
            '27-Nov','2-Dec','7-Dec',...
            '27-Nov','2-Dec','7-Dec',...
            '27-Nov','2-Dec','7-Dec',...
            '27-Nov','2-Dec','7-Dec',...
            '27-Nov','2-Dec','7-Dec',...
            '27-Nov','2-Dec','7-Dec',...
            '27-Nov','7-Dec',...
            });
end
xtickangle(45)
grid on
xlabel('sample')
ylabel('score')

figure
plot(mzs0(sum(isnan(data4pca))==0),coeff(:,[1 2]))

%%

[ ~, pc2_peaksOfInteresti ] = sort(abs(coeff(:,2)),'descend');

%%

figure;
hold on
data4pca = squeeze(mean(dataMatrix,1,'omitnan'))';
data4pca = (data4pca-mean(data4pca,'omitnan'))./std(data4pca(sum(~isnan(data4pca)==0,2)==0,:));
plot(mean(data4pca,2)', 'k', 'linewidth',3)
%plot(mean(squeeze(sum(dataMatrix,2)),1), 'b', 'linewidth',3)
xticks(1:size(dataMatrix,3))
switch data_type
    case 'vbh & tumours'
        xticklabels({...
            '27-Nov','2-Dec','7-Dec',...
            '27-Nov','2-Dec','7-Dec',...
            '27-Nov','2-Dec','7-Dec',...
            '27-Nov','2-Dec','7-Dec',...
            '27-Nov','2-Dec','7-Dec',...
            '27-Nov','2-Dec','7-Dec',...
            '27-Nov','7-Dec',...
            });
end
xtickangle(45)
grid on
xlabel('sample')
ylabel('TIC')