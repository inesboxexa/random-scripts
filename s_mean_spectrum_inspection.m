
files = dir('X:\Beatson\PI3K study (drug swiss)\DPO\uMALDI Neg\spectra details\');

for i = 3:size(files,1)
    
    fileName = files(i).name;
    
    if i==3
        load(['X:\Beatson\PI3K study (drug swiss)\DPO\uMALDI Neg\spectra details\' fileName '\tissue only\totalSpectrum_mzvalues.mat'])
        mzs = totalSpectrum_mzvalues';
        mean_spectra = zeros(size(totalSpectrum_mzvalues,2),size(files,1));
        total_spectra = mean_spectra;
        
        load(['X:\Beatson\PI3K study (drug swiss)\DPO\uMALDI Neg\spectra details\' fileName '\tissue only\peakDetails.mat'])
        [ ~, peaksi ] = sort(peakDetails(:,4),'descend');
        N = 2000;
        top_peaks_mzs = peakDetails(peaksi(1:N),2);
    
    end
    
    load(['X:\Beatson\PI3K study (drug swiss)\DPO\uMALDI Neg\spectra details\' fileName '\tissue only\totalSpectrum_intensities.mat'])
    load(['X:\Beatson\PI3K study (drug swiss)\DPO\uMALDI Neg\spectra details\' fileName '\tissue only\pixels_num.mat'])
        
    mean_spectra(:,i) = totalSpectrum_intensities'./pixels_num;
    total_spectra(:,i) = totalSpectrum_intensities';

end

mzsOfinterest = (mzs>=50).*(mzs<=1000);

mzs(~mzsOfinterest)=[];
mean_spectra(~mzsOfinterest,:)=[];
total_spectra(~mzsOfinterest,:)=[];

%%

figure; plot(mzs,mean_spectra); hold on; stem(top_peaks_mzs,repmat(max(mean_spectra(:)),size(top_peaks_mzs,1),1),'k')