
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

%% Gathering the veal-brain-homogenate data

% Please list the paths to all folders containing data.

data_folders = { 'D:\veal-brain-homogenate-study\data\' };

% Please list the strings that matches the names of the files to be analised. Each row should match each folder specified above. If all files need to be analised, please use '*'.

dataset_name_portion = { '*' };

% ! Please don't modify the code from here till end of this cell.

filesToProcess = []; for i = 1:length(data_folders); filesToProcess = [ filesToProcess; dir([data_folders{i} filesep dataset_name_portion{i} '.imzML']) ]; end % Files and adducts information gathering

% Defining the group of files of interest

MyFilesInfo = filesToProcess([ 13 16 17 21 33 38 ]);

dataMatrix = [];
sample_ids = [];
figure;
for i = 1:size(MyFilesInfo,1)
    
    file_name = strrep(MyFilesInfo(i).name(1:end-6), '_','-');
    
    load([ 'D:\veal-brain-homogenate-study\dpo\spectra details' filesep MyFilesInfo(i).name(1:end-6) '\no mask\' char(normalisation) filesep 'data' ]);
    load([ 'D:\veal-brain-homogenate-study\dpo\spectra details' filesep MyFilesInfo(i).name(1:end-6) '\no mask\datacubeonly_peakDetails' ]);
    load([ 'D:\veal-brain-homogenate-study\dpo\spectra details' filesep MyFilesInfo(i).name(1:end-6) '\no mask\width' ]);
    load([ 'D:\veal-brain-homogenate-study\dpo\spectra details' filesep MyFilesInfo(i).name(1:end-6) '\no mask\height' ]);
    load([ 'D:\veal-brain-homogenate-study\dpo\spectra details' filesep MyFilesInfo(i).name(1:end-6) '\no mask\totalSpectrum_intensities' ]);
    load([ 'D:\veal-brain-homogenate-study\dpo\spectra details' filesep MyFilesInfo(i).name(1:end-6) '\no mask\totalSpectrum_mzvalues' ]);
    
    dataMatrix0 = data;
    mzs0 = datacubeonly_peakDetails(:,2);
    
    if isempty(dataMatrix)
        dataMatrix = dataMatrix0;
        mzs = mzs0;
        spectraMatrix = totalSpectrum_intensities';
        spectra_mzs = totalSpectrum_mzvalues';
    else
        dataMatrix = cat(3,dataMatrix,dataMatrix0);
        mzs = [mzs,mzs0];
        spectraMatrix = [spectraMatrix,totalSpectrum_intensities'];
    end
                        
    % cv = std(dataMatrix0,1)./mean(dataMatrix0,1);
    subplot(2,4,i)
    tic_image = reshape(sum(data,2),width,height)';
    imagesc(tic_image); axis image off; colorbar
    title(file_name)
    
    sample_ids = [ sample_ids, string(file_name) ];
            
end

%%

figure
plot(spectra_mzs,spectraMatrix);
legend({ '20201126-DESI-NEG-VBH-2x2-3','20201127-DESI-NEG-VBH-2x2-3','20201130-DESI-NEG-VBH-2x2','20201202-DESI-NEG-VBH-2x2-3','20201207-xDESI-V3-NEG-VBH','20201209-DESI-NEG-VBH-2x2-5' })

%% Gathering the tumours data for the same dates

data_type = 'tumours';

% Please list the paths to all folders containing data.

data_folders = { 
    'D:\veal-brain-homogenate-study\tumours data\',...
    };

% Please list the strings that matches the names of the files to be analised. Each row should match each folder specified above. If all files need to be analised, please use '*'.

dataset_name_portion = { '*' };

% ! Please don't modify the code from here till end of this cell.

filesToProcess = []; for i = 1:length(data_folders); filesToProcess = [ filesToProcess; dir([data_folders{i} filesep dataset_name_portion{i} '.imzML']) ]; end % Files and adducts information gathering

% Defining the group of files of interest

MyFilesInfo = filesToProcess([ ... ]);

dataMatrix = [];
sample_ids = [];
figure;
for i = 1:size(MyFilesInfo,1)
    
    file_name = strrep(MyFilesInfo(i).name(1:end-6), '_','-');
    
    load([ 'D:\veal-brain-homogenate-study\tumours dpo\spectra details' filesep MyFilesInfo(i).name(1:end-6) '\tissue only\' char(normalisation) filesep 'data' ]);
    load([ 'D:\veal-brain-homogenate-study\tumours dpo\spectra details' filesep MyFilesInfo(i).name(1:end-6) '\tissue only\datacubeonly_peakDetails' ]);
    load([ 'D:\veal-brain-homogenate-study\tumours dpo\spectra details' filesep MyFilesInfo(i).name(1:end-6) '\tissue only\width' ]);
    load([ 'D:\veal-brain-homogenate-study\tumours dpo\spectra details' filesep MyFilesInfo(i).name(1:end-6) '\tissue only\height' ]);
    load([ 'D:\veal-brain-homogenate-study\tumours dpo\spectra details' filesep MyFilesInfo(i).name(1:end-6) '\tissue only\totalSpectrum_intensities' ]);
    load([ 'D:\veal-brain-homogenate-study\tumours dpo\spectra details' filesep MyFilesInfo(i).name(1:end-6) '\tissue only\totalSpectrum_mzvalues' ]);
    
    dataMatrix0 = data;
    mzs0 = datacubeonly_peakDetails(:,2);
    
    if isempty(dataMatrix)
        dataMatrix = dataMatrix0;
        mzs = mzs0;
        spectraMatrix = totalSpectrum_intensities';
        spectra_mzs = totalSpectrum_mzvalues';
    else
        if size(dataMatrix,1)<size(dataMatrix0,1)
            dataMatrix = cat(1,dataMatrix,NaN*ones(size(dataMatrix0,1)-size(dataMatrix,1),size(dataMatrix,2),size(ataMatrix,3)));
        elseif size(dataMatrix,1)>size(dataMatrix0,1)
            dataMatrix0 = cat(1,dataMatrix0,NaN*ones(size(dataMatrix,1)-size(dataMatrix0,1),size(dataMatrix,2),size(ataMatrix,3))
        end
        
        dataMatrix = cat(3,dataMatrix,dataMatrix0);
        mzs = [mzs,mzs0];
        spectraMatrix = [spectraMatrix,totalSpectrum_intensities'];
    end
                        
    % cv = std(dataMatrix0,1)./mean(dataMatrix0,1);
    subplot(2,4,i)
    tic_image = reshape(sum(data,2),width,height)';
    imagesc(tic_image); axis image off; colorbar
    title(file_name)
    
    sample_ids = [ sample_ids, string(file_name) ];
            
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

peaksOfInterest_lowestp = [
    500.278
    235.168
    244.136
    290.1
    736.644
    679.504
    737.646
    460.304
    680.508
    582.28
    329.192
    677.49
    690.532
    528.308
    356.262
    897.752
    717.562
    731.114
    440.264
    973.782
    734.56
    708.61
    ];

peaksOfInterest_largestp = [
   168.99
   296.014
   442.016
   462.964
   465.988
   481.964
   501.946
   504.346
   507.992
   517.96
   660.256
   958.308
    ];

%%

peaksOfInterest = peaksOfInterest_largestp;

peaksOfInteresti = [];
for peak = peaksOfInterest'
    [~, index] = min(abs(mzs0-peak));
    peaksOfInteresti = [ peaksOfInteresti; index];
end

% peaksOfInteresti = pc2_peaksOfInteresti(1:12);

% Update data and spectral channels

dataMatrix2plot = dataMatrix(:,peaksOfInteresti,:);
mzs2plot = mzs0(peaksOfInteresti);

% % Remove peaks with zero meadian across data
%
% [ peaks2keep, peaks2keepi ] = sort(prod(squeeze(median(dataMatrix,1,'omitnan')),2),'descend');
% dataMatrix = dataMatrix(:,peaks2keepi(peaks2keep>0),:);
%
% mzs = mzs(peaks2keepi(peaks2keep>0),1);
%
% % Select random peaks
%
% % randpeaksi = sort(randi(size(dataMatrix,2),[1 3*3]));
% % randpeaksmz = mzs(randpeaksi,1);
%
% randpeaksi = 1:12;
% randpeaksmz = mzs(randpeaksi,1);

means = squeeze(mean(dataMatrix2plot(:,:,:),1,'omitnan'));
medians = squeeze(median(dataMatrix2plot(:,:,:),1,'omitnan'));
prctile75 = squeeze(prctile(dataMatrix2plot(:,:,:),75,1));
prctile25 = squeeze(prctile(dataMatrix2plot(:,:,:),25,1));

colors = viridis(size(means,1));

figure
for peaki = 1:size(means,1)
    subplot(3,4,peaki)
    hold on
    title(['m/z ',  num2str(round(mzs2plot(peaki),3))])
    % plot(1:size(totalcount,2),totalcount(peaki,:),'-o','color', colors(peaki,:),'MarkerEdgeColor',colors(peaki,:),'MarkerFaceColor',colors(peaki,:))
    plot(1:size(means,2),means(peaki,:),'-x','color', colors(peaki,:),'MarkerEdgeColor',colors(peaki,:),'MarkerFaceColor',[1 1 1])
    plot(1:size(means,2),medians(peaki,:),'-o','color', colors(peaki,:),'MarkerEdgeColor',colors(peaki,:),'MarkerFaceColor',colors(peaki,:))
    plot(1:size(means,2),prctile75(peaki,:),':v','color', colors(peaki,:),'MarkerEdgeColor',colors(peaki,:),'MarkerFaceColor',[1 1 1])
    plot(1:size(means,2),prctile25(peaki,:),':^','color', colors(peaki,:),'MarkerEdgeColor',colors(peaki,:),'MarkerFaceColor',[1 1 1])
    grid on
    ylabel('ion counts')
    xlabel('sample')
    xticks(1:size(means,2))
    switch data_type
        case 'vbh'
            xticklabels({'20201126','20201127','20201130','20201202','20201207','20201209'});
        case 'tumours'
            xticklabels({...
                '20201126','20201126','20201126','20201126','20201126','20201126','20201126',...
                '20201127','20201127','20201127','20201127',...
                '20201130',...
                '20201202','20201202','20201202','20201202',...
                '20201207','20201207','20201207',...
                '20201209','20201209','20201209','20201209','20201209'
                });
    end
    xtickangle(45)
    if peaki == 4 % size(totalcount,1)
    legend({'mean','median','percentile 75','percentile 25'})
    legend('boxoff')
    legend('Location','best')
    end
end

%% Running a Principal Component Analysis (PCA) of the 500 most intense ions and plotting the scores which illustrate the intensaty trends across samples

colors = tab10;
figure
data4pca = squeeze(mean(dataMatrix,1,'omitnan'))';
data4pca = (data4pca-mean(data4pca,'omitnan'))./std(data4pca(sum(~isnan(data4pca)==0,2)==0,:));
[coeff,score,latent,tsquared,explained,mu] = pca(data4pca);
bar(score(:,1:5),'EdgeColor',[1 1 1])
% stem(1:size(data4pca,1),score(:,1),'-o','color', colors(1,:),'MarkerEdgeColor',colors(1,:),'MarkerFaceColor',colors(1,:))
% stem(1:size(data4pca,1),score(:,2),'-o','color', colors(3,:),'MarkerEdgeColor',colors(3,:),'MarkerFaceColor',colors(3,:))
% stem(1:size(data4pca,1),score(:,3),'-o','color', colors(5,:),'MarkerEdgeColor',colors(5,:),'MarkerFaceColor',colors(5,:))
% stem(1:size(data4pca,1),score(:,4),'-o','color', colors(7,:),'MarkerEdgeColor',colors(7,:),'MarkerFaceColor',colors(7,:))
% stem(1:size(data4pca,1),score(:,5),'-o','color', colors(9,:),'MarkerEdgeColor',colors(9,:),'MarkerFaceColor',colors(9,:))
legend({['PC 1 (', num2str(round(explained(1),2)) '% var expl)'], ['PC 2 (', num2str(round(explained(2),2)) '% var expl)'], ['PC 3 (', num2str(round(explained(3),2)) '% var expl)'], ['PC 4 (', num2str(round(explained(4),2)) '% var expl)'], ['PC 5 (', num2str(round(explained(5),2)) '% var expl)']});
legend('boxoff')
legend('Location','best')
xticks(1:size(means,2))
xticklabels({'20201126','20201127','20201130','20201202','20201207','20201209'});
xtickangle(45)
grid on
xlabel('sample')
ylabel('score')

figure
plot(mzs0,coeff(:,[1 2]))

%%

[ ~, pc2_peaksOfInteresti ] = sort(abs(coeff(:,2)),'descend');
