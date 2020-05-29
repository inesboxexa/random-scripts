
threshold = 100; % threshold for the mean intensity

peak_assigments_folder = 'X:\PDAC Combo\az neg desi\dpo\peak assignments\';

filefolder_struct = dir(peak_assigments_folder);

for i = 3:size(filefolder_struct)
    
    filefolder_name = filefolder_struct(i).name;
    
    path = [peak_assigments_folder filesep filefolder_name '\tissue only\']; % location of the assigments mat file
    file = 'relevant_lists_sample_info'; % name of the assigments mat file
    
    load([ path filesep file ])
    
    % save([path filesep 'original_' file],'relevant_lists_sample_info','-v7.3') % save original assignments mat file
    
    mean_intensity = double(relevant_lists_sample_info(:,11)); % collect mean intensities array
    
    relevant_lists_sample_info(mean_intensity<threshold,:) = []; % remove mean intensities below threshold
    
    save([path filesep file],'relevant_lists_sample_info','-v7.3') % save new assignments mat file
    
end