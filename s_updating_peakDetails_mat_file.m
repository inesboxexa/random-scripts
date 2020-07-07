
dpo_path = 'X:\Teresa\pancreatic tumour - blue slices set\dpo\';
mask = 'no mask';

spectra_details_folder = [ dpo_path 'spectra details\'];
peak_assigments_folder = [ dpo_path 'peak assignments\'];

filefolder_struct = dir(peak_assigments_folder);

for i = 3:size(filefolder_struct)
    
    filefolder_name = filefolder_struct(i).name;
    
    % Total Spectrum Intensities
    
    path = [spectra_details_folder filesep filefolder_name '\' mask '\']; % location of the total spectrum mat file
    file = 'totalSpectrum_intensities'; % name of the total spectrum mat file
    
    load([ path filesep file ])
    
    % threshold = 100; % threshold for the mean intensity
    
    % 3 times the median of the total intensities threshold
    
    threshold = 3*median(totalSpectrum_intensities(totalSpectrum_intensities>0)); % 3 times the median of the total intensities above 0
    
    % Peak Details
    
    path = [spectra_details_folder filesep filefolder_name '\' mask '\']; % location of the assigments mat file
    file = 'peakDetails'; % name of the assigments mat file
    
    load([ path filesep file ])
    
    figure('Name',filefolder_name); 
    plot(totalSpectrum_mzvalues, totalSpectrum_intensities, 'k');
    hold on; stem(totalSpectrum_mzvalues(totalSpectrum_intensities>threshold),totalSpectrum_intensities(totalSpectrum_intensities >= threshold),'.r')
    
    save([path filesep 'original_' file],'relevant_lists_sample_info','-v7.3') % save original assignments mat file
    
    % Mean intensity based threshold
    
    % mean_intensity = double(relevant_lists_sample_info(:,11)); % mean intensities array
    % relevant_lists_sample_info(mean_intensity<threshold,:) = []; % remove mean intensities below threshold
    
    % 3x the median of the total intensities threshold
    
    total_intensity = double(peakDetails(:,4)); % total intensities array
    relevant_lists_sample_info(total_intensity < threshold,:) = []; % remove mean intensities below threshold
    
    save([path filesep file],'relevant_lists_sample_info','-v7.3') % save new assignments mat file
    
end