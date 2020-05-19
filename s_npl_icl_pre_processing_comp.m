
% cd('X:\Alex\Data study\')

load('X:\Alex\Data study\ICL datacubes\ICL-NoMask.mat')
load('X:\Alex\Data study\20181029_PDAC_Combo_set1_slide11_A1_desi_neg.mat')
load('X:\Alex\Data study\peakMathcingResults.mat')

npl_data = data;
npl_mzvalues = peaks1';
npl_mzindexes = {}; for i = 1:length(commonIdx1); npl_mzindexes{i} = commonIdx1{i}'; end

icl_data = icl.sp';
icl_mzvalues = peaks2;
icl_mzindexes = {}; for i = 1:length(commonIdx2); icl_mzindexes{i} = commonIdx2{i}'; end

width = icl.sz(1);
height = icl.sz(2);

%% Single Ion Images Comparison

ppm_list = [5, 10, 25, 50, 100];

for set = 1:length(npl_mzindexes)
    
    random_mz4fig = randperm(length(npl_mzindexes{set}),100);
    
    images_sim_metrics = NaN*ones(6,length(npl_mzindexes{set}));
    scaled_images_sim_metrics = NaN*ones(6,length(npl_mzindexes{set}));
    
    mkdir(['X:\Alex\Data study\NPL to ICL comparison\' num2str(ppm_list(set)) 'ppm tolerance\image similarity metrics\'])
    cd(['X:\Alex\Data study\NPL to ICL comparison\' num2str(ppm_list(set)) 'ppm tolerance\image similarity metrics\'])
    
    for mz_index = 1:length(npl_mzindexes{set})
        
        npl_ion_image = npl_data(:,npl_mzindexes{set}(mz_index));
        icl_ion_image = icl_data(:,icl_mzindexes{set}(mz_index));
        
        scaled_npl_ion_image = (npl_ion_image-min(npl_ion_image(:)))./max(npl_ion_image(:)-min(npl_ion_image(:)));
        scaled_icl_ion_image = (icl_ion_image-min(icl_ion_image(:)))./max(icl_ion_image(:)-min(icl_ion_image(:)));
        
        % Correlation
        
        images_sim_metrics(1, mz_index) = corr(npl_ion_image,icl_ion_image,'type','Pearson');
        images_sim_metrics(2, mz_index) = corr(npl_ion_image,icl_ion_image,'type','Kendall');
        images_sim_metrics(3, mz_index) = corr(npl_ion_image,icl_ion_image,'type','Spearman');
        
        scaled_images_sim_metrics(1, mz_index) = corr(scaled_npl_ion_image,scaled_icl_ion_image,'type','Pearson');
        scaled_images_sim_metrics(2, mz_index) = corr(scaled_npl_ion_image,scaled_icl_ion_image,'type','Kendall');
        scaled_images_sim_metrics(3, mz_index) = corr(scaled_npl_ion_image,scaled_icl_ion_image,'type','Spearman');
        
        % Structural Similarity Index (SSIM)
        
        images_sim_metrics(4, mz_index) = ssim(npl_ion_image,icl_ion_image);
        scaled_images_sim_metrics(4, mz_index) = ssim(scaled_npl_ion_image,scaled_icl_ion_image);
        
        % Peak signal-to-noise ratio (pSNR)
        
        images_sim_metrics(5, mz_index) = psnr(npl_ion_image,icl_ion_image);
        scaled_images_sim_metrics(5, mz_index) = psnr(scaled_npl_ion_image,scaled_icl_ion_image);
        
        % Mean-squared error (MSE)
        
        images_sim_metrics(6, mz_index) = immse(npl_ion_image,icl_ion_image);
        scaled_images_sim_metrics(6, mz_index) = immse(scaled_npl_ion_image,scaled_icl_ion_image);
        
        if sum(mz_index == random_mz4fig) == 1
            
            fig = figure('units','normalized','outerposition',[0 0 1 1]);
            
            subplot(2,2,1)
            imagesc(reshape(npl_ion_image,width,height)'); axis image off; colormap viridis; colorbar
            title({['NPL - mz ' num2str(npl_mzvalues(npl_mzindexes{set}(mz_index)))]},'FontSize',16)
            
            subplot(2,2,3)
            imagesc(reshape(icl_ion_image,width,height)'); axis image off; colormap viridis; colorbar
            title({['ICL - mz ' num2str(icl_mzvalues(icl_mzindexes{set}(mz_index)))]},'FontSize',16)
            
            subplot(2,2,2)
            I = ones(height,width,3);
            I(:,:,1) = 1-reshape(scaled_npl_ion_image,width,height)';
            I(:,:,2) = 1-reshape(scaled_icl_ion_image,width,height)';
            mask = repmat(logical(sum(I,3)==0),1,1,3);
            I(mask) = 1;
            imagesc(I); axis image off;
            title({'NPL (cian) - ICL (magenta) (intenisties are scaled between 0 and 1)'},'FontSize',16)
            
            subplot(2,2,4)
            cdata = [
                images_sim_metrics(:, mz_index)'
                scaled_images_sim_metrics(:, mz_index)'
                ];
            xvalues = {'Pearson corr', 'Kendall corr', 'Spearman corr', 'ssim', 'pSNR', 'MSE'};
            yvalues = {'Original','Scaled'};
            h = heatmap(xvalues,yvalues,cdata);
            h.Title = 'Similarity Metrics';
            h.XLabel = 'Metric';
            h.YLabel = 'Data';
            h.ColorScaling = 'scaledcolumns';
            
            savefig(fig,['image_overlay_NPL_mz' num2str(npl_mzvalues(npl_mzindexes{set}(mz_index))) '_ICL_mz' num2str(icl_mzvalues(icl_mzindexes{set}(mz_index))) '.fig'])
            saveas(fig,['image_overlay_NPL_mz' num2str(npl_mzvalues(npl_mzindexes{set}(mz_index))) '_ICL_mz' num2str(icl_mzvalues(icl_mzindexes{set}(mz_index))) '.png'])
            
            close all
            
        end
        
    end
    
    save('images_sim_metrics','images_sim_metrics')
    save('scaled_images_sim_metrics','scaled_images_sim_metrics')
    
end

%% Metrics by mass

metric = {'Pearson corr', 'Kendall corr', 'Spearman corr', 'ssim', 'pSNR', 'MSE'};

ppm_list = [ 5, 10, 25, 50, 100 ];

for set = 1:length(npl_mzindexes)
    
    cd(['X:\Alex\Data study\NPL to ICL comparison\' num2str(ppm_list(set)) 'ppm tolerance\image similarity metrics\'])
    
    load('images_sim_metrics')
    load('scaled_images_sim_metrics')
    
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2,3,1:3)
    plot(npl_mzvalues(npl_mzindexes{set}),images_sim_metrics(2,:)','o','color',[.85 .85 .85]); hold on
    plot(npl_mzvalues(npl_mzindexes{set}),images_sim_metrics(3,:)','d','color',[.85 .85 .85]); hold on
    plot(npl_mzvalues(npl_mzindexes{set}),images_sim_metrics(1,:)','xk'); hold on
    legend({'Kendall corr', 'Spearman corr','Pearson corr'})
    grid on
    for metric_i = 4:6
        subplot(2,3,metric_i)
        plot(npl_mzvalues(npl_mzindexes{set}),[images_sim_metrics(metric_i,:)' -scaled_images_sim_metrics(metric_i,:)'],'.')
        title({metric{metric_i}})
        legend({'Original','Scaled'})
        grid on
    end
    
    savefig(fig,'metrics_vs_mass.fig')
    
    close all
    
end

%% Metrics by mass

ppm_list = [ 5, 10, 25, 50, 100 ];

for set = 1:length(npl_mzindexes)
    
    cd(['X:\Alex\Data study\NPL to ICL comparison\' num2str(ppm_list(set)) 'ppm tolerance\npl-icl mz differences\'])
    
    for mz_index = 1:length(npl_mzindexes{set})
        
        npl_ion_image = npl_data(:,npl_mzindexes{set}(mz_index));
        icl_ion_image = icl_data(:,icl_mzindexes{set}(mz_index));
        
        
        fig = figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(2,3,1:3)
        plot(npl_mzvalues(npl_mzindexes{set}),images_sim_metrics(2,:)','o','color',[.85 .85 .85]); hold on
        plot(npl_mzvalues(npl_mzindexes{set}),images_sim_metrics(3,:)','d','color',[.85 .85 .85]); hold on
        plot(npl_mzvalues(npl_mzindexes{set}),images_sim_metrics(1,:)','xk'); hold on
        legend({'Kendall corr', 'Spearman corr','Pearson corr'})
        grid on
        for metric_i = 4:6
            subplot(2,3,metric_i)
            plot(npl_mzvalues(npl_mzindexes{set}),[images_sim_metrics(metric_i,:)' -scaled_images_sim_metrics(metric_i,:)'],'.')
            title({metric{metric_i}})
            legend({'Original','Scaled'})
            grid on
        end
        
        savefig(fig,'metrics_vs_mass.fig')
        
        close all
        
    end
    
end

%% Drugs and related metabolites

% Brain

load('C:\Users\tm6\Documents\HCP Anywhere\CRUK Data processing study\CassetteDosed\npl data processing outputs (mat files)\neg_brain_data.mat')
load('C:\Users\tm6\Documents\HCP Anywhere\CRUK Data processing study\CassetteDosed\npl data processing outputs (mat files)\neg_brain_mzvalues.mat')
load('C:\Users\tm6\Documents\HCP Anywhere\CRUK Data processing study\CassetteDosed\npl data processing outputs (mat files)\neg_brain_width.mat')
load('C:\Users\tm6\Documents\HCP Anywhere\CRUK Data processing study\CassetteDosed\npl data processing outputs (mat files)\neg_brain_height.mat')

npl_data = neg_brain_data;
npl_mzvalues = neg_brain_mzvalues;
npl_width = neg_brain_width;
npl_height = neg_brain_height;


load('C:\Users\tm6\Documents\HCP Anywhere\CRUK Data processing study\CassetteDosed\icl data processing outputs (h5 and mat files)\ICL-Brain-Neg.mat')

icl_data = icl.sp;
icl_mzvalues = icl.mz;
icl_width = icl.sz(1);
icl_height = icl.sz(2);

mkdir('D:\icl-npl data study\CassetteDosed\neg-brain-sii\')
cd('D:\icl-npl data study\CassetteDosed\neg-brain-sii\')

%%

% Kidney

load('C:\Users\tm6\Documents\HCP Anywhere\CRUK Data processing study\CassetteDosed\npl data processing outputs (mat files)\neg_kidney_data.mat')
load('C:\Users\tm6\Documents\HCP Anywhere\CRUK Data processing study\CassetteDosed\npl data processing outputs (mat files)\neg_kidney_mzvalues.mat')
load('C:\Users\tm6\Documents\HCP Anywhere\CRUK Data processing study\CassetteDosed\npl data processing outputs (mat files)\neg_kidney_width.mat')
load('C:\Users\tm6\Documents\HCP Anywhere\CRUK Data processing study\CassetteDosed\npl data processing outputs (mat files)\neg_kidney_height.mat')

npl_data = neg_kidney_data;
npl_mzvalues = neg_kidney_mzvalues;
npl_width = npl_kidney_width;
npl_height = npl_kidney_height;


load('C:\Users\tm6\Documents\HCP Anywhere\CRUK Data processing study\CassetteDosed\icl data processing outputs (h5 and mat files)\ICL-Kidney-Neg.mat')

icl_data = icl.sp;
icl_mzvalues = icl.mz;
icl_width = icl.sz(1);
icl_height = icl.sz(2);

mkdir('D:\icl-npl data study\CassetteDosed\neg-kidney-sii\')
cd('D:\icl-npl data study\CassetteDosed\neg-kidney-sii\')

%% PDAC

% drugs = [ 475.0927, 454.9927, 493.1033, 489.1823, 461.2307, 479.2412, 453.2056, 497.2073, 443.2201, 435.195, 471.2162, 457.0821, 511.0694, 436.9822, 436.9822, 473.0033 ];

% Dossed Cassette

% % pos adducts
% 
% drugs = [ ...
%     401.175085, 402.182361, 424.164306, 440.138244, 384.171796, ... % Moxifloxacin
%     312.140868, 313.148144, 335.130089, 351.104027, 295.137579, ... % Olanzapine
%     393.168857, 394.176133, 416.158078, 432.132016, 376.165568, ... % Erlotinib
%     471.313729, 472.321005, 494.30295, 510.276888, 454.31044 ... % Terfenadine
%     ];
% 
% mkdir('X:\Alex\Data study\CassetteDosed\drugs\pos\')
% cd('X:\Alex\Data study\CassetteDosed\drugs\pos\')

%%

% neg adducts

drugs = [ ...
    392.1616, 500.2806, 410.1721, 396.1565, 378.1459, ...
    506.2831, 428.1383, 436.1445, 333.0946, 518.2912, ...
    414.1226, 482.2701, 311.1336, 297.1179, 347.1103, ...
    315.1285, 315.1285, 418.1784, 400.1678, 418.1784, ...
    400.1678, 400.1678, 488.317, 360.1354, 360.1354, ...
    504.6761, 486.6656, 522.6422, 522.6422, 468.655, ...
    468.655, 363.1051 ...
    ];

%% Single ion images for drugs!

for mz = drugs(end)
    
    [ ~, npl_mz_index ] = min(abs(npl_mzvalues-mz));
    npl_ion_image = npl_data(:,npl_mz_index);
    
    [ ~, icl_mz_index ] = min(abs(icl_mzvalues-mz));
    icl_ion_image = icl_data(:,icl_mz_index);
    
    scaled_npl_ion_image = (npl_ion_image-min(npl_ion_image(:)))./max(npl_ion_image(:)-min(npl_ion_image(:)));
    scaled_icl_ion_image = (icl_ion_image-min(icl_ion_image(:)))./max(icl_ion_image(:)-min(icl_ion_image(:)));
    
%     % Correction the different sizes!!!!
%     
%     icl_ion_image = icl_ion_image(1:end-icl_width,:);
%     scaled_icl_ion_image = scaled_icl_ion_image(1:end-icl_width,:);
    
    width = npl_width;
    height = npl_height;
    
    % Correlation
    
    images_sim_metrics(1, 1) = corr(npl_ion_image,icl_ion_image,'type','Pearson');
    images_sim_metrics(2, 1) = corr(npl_ion_image,icl_ion_image,'type','Kendall');
    images_sim_metrics(3, 1) = corr(npl_ion_image,icl_ion_image,'type','Spearman');
    
    scaled_images_sim_metrics(1, 1) = corr(scaled_npl_ion_image,scaled_icl_ion_image,'type','Pearson');
    scaled_images_sim_metrics(2, 1) = corr(scaled_npl_ion_image,scaled_icl_ion_image,'type','Kendall');
    scaled_images_sim_metrics(3, 1) = corr(scaled_npl_ion_image,scaled_icl_ion_image,'type','Spearman');
    
    % Structural Similarity Index (SSIM)
    
    images_sim_metrics(4, 1) = ssim(npl_ion_image,icl_ion_image);
    scaled_images_sim_metrics(4, 1) = ssim(scaled_npl_ion_image,scaled_icl_ion_image);
    
    % Peak signal-to-noise ratio (pSNR)
    
    images_sim_metrics(5, 1) = psnr(npl_ion_image,icl_ion_image);
    scaled_images_sim_metrics(5, 1) = psnr(scaled_npl_ion_image,scaled_icl_ion_image);
    
    % Mean-squared error (MSE)
    
    images_sim_metrics(6, 1) = immse(npl_ion_image,icl_ion_image);
    scaled_images_sim_metrics(6, 1) = immse(scaled_npl_ion_image,scaled_icl_ion_image);
    
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    
    subplot(2,2,1)
    imagesc(reshape(npl_ion_image,width,height)'); axis image off; colormap viridis; colorbar
    title({['NPL - mz ' num2str(npl_mzvalues(npl_mz_index)) ' / npl-theo error ' num2str((mz-npl_mzvalues(npl_mz_index))./(npl_mzvalues(npl_mz_index)).*1000000) ' / npl-icl error ' num2str((npl_mzvalues(npl_mz_index)-icl_mzvalues(icl_mz_index))./(npl_mzvalues(npl_mz_index)).*1000000)]},'FontSize',14)
    
    subplot(2,2,3)
    imagesc(reshape(icl_ion_image,width,height)'); axis image off; colormap viridis; colorbar
    title({['ICL - mz ' num2str(icl_mzvalues(icl_mz_index)) ' / icl-theo error ' num2str((mz-icl_mzvalues(icl_mz_index))./(icl_mzvalues(icl_mz_index)).*1000000) ' / npl-icl error ' num2str((npl_mzvalues(npl_mz_index)-icl_mzvalues(icl_mz_index))./(npl_mzvalues(npl_mz_index)).*1000000)]},'FontSize',14)
    
    subplot(2,2,2)
    I = ones(height,width,3);
    I(:,:,1) = 1-reshape(scaled_npl_ion_image,width,height)';
    I(:,:,2) = 1-reshape(scaled_icl_ion_image,width,height)';
    mask = repmat(logical(sum(I,3)==0),1,1,3);
    I(mask) = 1;
    imagesc(I); axis image off;
    title({'NPL (cian) - ICL (magenta) (intenisties are scaled between 0 and 1)'},'FontSize',16)
    
    subplot(2,2,4)
    cdata = [
        images_sim_metrics'
        scaled_images_sim_metrics'
        ];
    xvalues = {'Pearson corr', 'Kendall corr', 'Spearman corr', 'ssim', 'pSNR', 'MSE'};
    yvalues = {'Original','Scaled'};
    h = heatmap(xvalues,yvalues,cdata);
    h.Title = 'Similarity Metrics';
    h.XLabel = 'Metric';
    h.YLabel = 'Data';
    h.ColorScaling = 'scaledcolumns';
    
    savefig(fig,['image_overlay_NPL_mz' num2str(npl_mzvalues(npl_mz_index)) '_ICL_mz' num2str(icl_mzvalues(icl_mz_index)) '.fig'])
    saveas(fig,['image_overlay_NPL_mz' num2str(npl_mzvalues(npl_mz_index)) '_ICL_mz' num2str(icl_mzvalues(icl_mz_index)) '.png'])
    
    close all
    
end