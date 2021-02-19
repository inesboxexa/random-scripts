cd('X:\Beatson\Intracolonic tumour study\dpo\neg DESI\rois\')

filesi = dir('*_R1B1_*');

for i = 1:size(filesi,1)
    
    cd([filesi(i).folder filesep filesi(i).name])
    
    filesii = dir();
    
    for ii = 3:size(filesii,1)
        
        cd([filesii(ii).folder filesep filesii(ii).name])
        
        load('roi')
        
        binary_mask = roi.pixelSelection;
        
        save('binary_mask','binary_mask')
        
    end
    
end