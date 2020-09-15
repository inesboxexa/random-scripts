
cd('X:\Beatson\PI3K study (drug swiss)\DPO\DESI Neg (Synapt)\whole tissue\rois\')

folders = dir; folders = folders(3:end,:);

for f = 1:size(folders,1)
   
    cd([ folders(f).folder filesep folders(f).name ])
    
    subfolders = dir; subfolders = subfolders(3:end,:);
    
    for sf = 1:size(subfolders,1)
        
        main_mask = 0;
        for clusterid = [ 1 2 5 6 ]
            
            cd([subfolders(sf).folder filesep 'mva based rois\PI3K SI DESI Negative\tissue only\tsne 7 components\RMS\mva 4000 highest peaks (2923 peaks discarded)\component ' num2str(clusterid)])
            load roi
            
            main_mask = main_mask + roi.pixelSelection;
            
        end
        
        roi = RegionOfInterest(roi.width,roi.height);
        roi.addPixels(main_mask)
        
        mkdir([subfolders(sf).folder filesep 'epithelium'])
        cd([subfolders(sf).folder filesep 'epithelium'])
        
        save('roi','roi')
        
        if sum(subfolders(sf).name == '_')>0
            
            clear roi
           
            cd([subfolders(sf).folder filesep subfolders(sf).name])
            load roi
            
            new_mask = main_mask.*roi.pixelSelection;
            
            roi = RegionOfInterest(roi.width,roi.height);
            roi.addPixels(new_mask)
            
            mkdir([subfolders(sf).folder filesep subfolders(sf).name ' epithelium'])
            cd([subfolders(sf).folder filesep subfolders(sf).name ' epithelium'])
            
            save('roi','roi')
                                    
        end
        
    end
    
end