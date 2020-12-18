
level1 = dir('X:\ICR Breast PDX\Data Processing Outputs\3d pdx pilot\rois\');

for i = 3:size(level1,1)    
    cd([ level1(i).folder filesep level1(i).name filesep 'mva based rois\neg desi 3d pdx pilot\tissue only\kmeans 8 components\no norm\mva 4000 highest peaks\'])
    level2 = dir('*component*');    
    for ii = 1:size(level2,1)        
        mkdir([level1(i).folder filesep level1(i).name filesep 'cluster-' level2(ii).name(11:end) '-' level1(i).name(18:21) '-' level1(i).name(23)])
        copyfile([level2(ii).folder filesep level2(ii).name filesep 'roi.mat'],[level1(i).folder filesep level1(i).name filesep 'cluster-' level2(ii).name(11:end) '-' level1(i).name(18:21) '-' level1(i).name(23)])
    end
end