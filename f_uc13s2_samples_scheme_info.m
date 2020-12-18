function [ extensive_filesToProcess, main_mask_list, smaller_masks_list, outputs_xy_pairs ] = f_uc13s2_samples_scheme_info( dataset_name, background, check_datacubes_size )

% output xy pairs determines position of each image, 1st is row 2nd is column
% to determine horizontal or vertical

switch dataset_name
    
    case "neg maldi UC13 S2 slide 1"
        
        data_folders = { 'X:\Hanifa\S2_GBM_August_2020_Data\MALDI\' };
        
        dataset_name = '*slide1*';
        
        filesToProcess = []; for i = 1:length(data_folders); filesToProcess = [ filesToProcess; dir([data_folders{i} dataset_name '.imzML']) ]; end
        
        if background == 1
            % with background
            main_mask_list = "no mask";
        else
            
            % tissue only
            
            main_mask_list = "tissue only";
            
            %
            
            extensive_filesToProcess(1:2,:) = filesToProcess(1:2,:);
            smaller_masks_list = [
                "Slide 1 brain coronal section";
                "Slide 1 brain oval section";
                ];
            
        end
        
        %
        
        outputs_xy_pairs = [
            1 1; 1 2;
            ];
        
     case "neg maldi UC13 S2 slide 1 oval"
        
        data_folders = { 'X:\Hanifa\S2_GBM_August_2020_Data\MALDI\' };
        
        dataset_name = '*Oval*';
        
        filesToProcess = []; for i = 1:length(data_folders); filesToProcess = [ filesToProcess; dir([data_folders{i} dataset_name '.imzML']) ]; end
        
        if background == 1
            % with background
            main_mask_list = "no mask";
        else
            
            % tissue only
            
            main_mask_list = "Slide 1 brain oval section";
            
            %
            
            extensive_filesToProcess(1,:) = filesToProcess(1,:);
            smaller_masks_list = [
                "Slide 1 brain oval section";
                ];
            
        end
        
        %
        
        outputs_xy_pairs = [
            1 1;
            ];
        
    case "neg maldi UC13 S2 slide 2"
        
        data_folders = { 'X:\Hanifa\S2_GBM_August_2020_Data\MALDI\' };
        
        dataset_name = '*slide2*';
        
        filesToProcess = []; for i = 1:length(data_folders); filesToProcess = [ filesToProcess; dir([data_folders{i} dataset_name '.imzML']) ]; end
        
        if background == 1
            % with background
            main_mask_list = "no mask";
        else
            
            % tissue only
            
            main_mask_list = "maldi slide2 oval";
            
            %
            
            extensive_filesToProcess(1,:) = filesToProcess(1,:);
            smaller_masks_list = [
                "maldi slide2 oval";
                ];
            
        end
        
        %
        
        outputs_xy_pairs = [
            1 1;
            ];
        
    case "neg desi UC13 S2 slide 1"
        
        data_folders = { 'X:\Hanifa\S2_GBM_August_2020_Data\DESI\' };
        
        dataset_name = '*slide1*';
        
        filesToProcess = []; for i = 1:length(data_folders); filesToProcess = [ filesToProcess; dir([data_folders{i} dataset_name '.imzML']) ]; end
        
        if background == 1
            % with background
            main_mask_list = "no mask";
        else
            
            % tissue only or ? other ROI
            
            main_mask_list = "oval rat brain";
            
            %
            
            extensive_filesToProcess(1,:) = filesToProcess(1,:);
            smaller_masks_list = [
                "oval rat brain";
                ];
            
        end
        
        %
        
        outputs_xy_pairs = [
            1 1;
            ];
        
    case "neg desi UC13 S2 slide 2"
        
        data_folders = { 'X:\Hanifa\S2_GBM_August_2020_Data\DESI\' };
        
        dataset_name = '*slide2*';
        
        filesToProcess = []; for i = 1:length(data_folders); filesToProcess = [ filesToProcess; dir([data_folders{i} dataset_name '.imzML']) ]; end
        
        if background == 1
            % with background
            main_mask_list = "no mask";
        else
            
            % tissue only
            
            main_mask_list = "desi slide2 oval";
            
            %
            
            extensive_filesToProcess(1,:) = filesToProcess(1,:);
            smaller_masks_list = [
                "desi slide2 oval";
                ];
            
        end
        
        %
        
        outputs_xy_pairs = [
            1 1;
            ];
        
        case "neg reims UC13 S2 slide 2"
        
        data_folders = { 'X:\Hanifa\S2_GBM_August_2020_Data\REIMS\Slide2\' };
        
        dataset_name = '*slide2*';
        
        filesToProcess = []; for i = 1:length(data_folders); filesToProcess = [ filesToProcess; dir([data_folders{i} dataset_name '.imzML']) ]; end
        
        if background == 1
            % with background
            main_mask_list = "no mask";
        else
            
            % tissue only
            
            main_mask_list = ["reims slide2 oval", "tissue only"];
            
            %
            
            extensive_filesToProcess(1,:) = filesToProcess(1,:);
            smaller_masks_list = [
                ["reims slide2 oval", "tissue only"];
                ];
            
        end
        
        %
        
        outputs_xy_pairs = [
            1 1;
            ];
        
        case "neg reims UC13 S2 slide 3 oval"
        
        data_folders = { 'X:\Hanifa\S2_GBM_August_2020_Data\REIMS\Slide3\' };
        
        dataset_name = '*oval*';
        
        filesToProcess = []; for i = 1:length(data_folders); filesToProcess = [ filesToProcess; dir([data_folders{i} dataset_name '.imzML']) ]; end
        
        if background == 1
            % with background
            main_mask_list = "no mask";
        else
            
            % tissue only
            
            main_mask_list = "reims slide3 oval";
            
            %
            
            extensive_filesToProcess(1,:) = filesToProcess(1,:);
            smaller_masks_list = [
                "reims slide3 oval";
                ];
            
        end
        
        %
        
        outputs_xy_pairs = [
            1 1;
            ];
        
        case "neg reims UC13 S2 slide 3 coronal"
        
        data_folders = { 'X:\Hanifa\S2_GBM_August_2020_Data\REIMS\Slide3\' };
        
        dataset_name = '*coronal*';
        
        filesToProcess = []; for i = 1:length(data_folders); filesToProcess = [ filesToProcess; dir([data_folders{i} dataset_name '.imzML']) ]; end
        
        if background == 1
            % with background
            main_mask_list = "no mask";
        else
            
            % tissue only
            
            main_mask_list = "reims slide3 coronal";
            
            %
            
            extensive_filesToProcess(1,:) = filesToProcess(1,:);
            smaller_masks_list = [
                "reims slide3 coronal";
                ];
            
        end
        
        %
        
        outputs_xy_pairs = [
            1 1;
            ];
end
