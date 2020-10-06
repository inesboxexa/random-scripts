% Script to remove peaks present in a spectrum A (with minimum predifined counts) from spectrum B, 
% written by Teresa Murta on 02 Oct 2020.

% Inputs: Path to one or more xlsx files with lists of mzs and counts for 2 spectrum, A and B.
% Outputs: One or more xlsx files with a clean version of spectrum B.

path = "X:\Teresa\background peak removal (Bins work)\"; % path to folder with xlsx files

files = [ "substrate peak remove data for TERESA.xlsx" ]; % names of all xlsx 

% Note: The xlsx file needs to have the background spectrum (mzs and counts) columns on
% the left side of the signals spectrum (mzs and counts).

for file = files
    
    [num,txt,raw] = xlsread([char(path) filesep char(file)]); % Load xlsx file
    
    num(:,isnan(sum(num,1))) = []; % remove columns of NaN
    
    background_mzs = num(:,1); % define background m/z values
    background_counts = num(:,2); % define background counts
    signal_mzs = num(:,3); % define signal m/z values
    signal_counts = num(:,4); % define signal counts

    clean_signal_counts = signal_counts; % define clean signal counts
    
    clean_background_mzs = background_mzs; % ignore background peaks with counts below 5e3
    clean_background_mzs(background_counts <= 5e3) = []; 
      
    peak_bol = 0;
    for background_peak = clean_background_mzs' % find peaks that are common to both signal and background
        [ minimum_ppm_error, signal_index ] = min(abs((signal_mzs-background_peak)./background_peak/1e6)); % minimum ppm error
        if minimum_ppm_error <= 3 % check if minimum ppm error is below 3
            clean_signal_counts(signal_index) = 0; % set signal counts to 0
        end
    end   
    
    char_file = char(file);
    xlswrite([char(path), filesep, char_file(1:end-5),' clean', char_file(end-4:end)], [signal_mzs clean_signal_counts]); % save clean signal peaks
    
    figure('Name', file); % plot
    stem(signal_mzs, signal_counts,'xk'); 
    hold on; stem(signal_mzs, clean_signal_counts,'xb'); 
    hold on; stem(background_mzs, -background_counts, 'xr');
    legend({'original signal', 'clean signal', 'background'});
    xlabel('m/z')
    ylabel('counts')
    grid on

end


