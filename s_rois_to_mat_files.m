
% read rois

path = 'X:\Alan\AnnotationTransfer\ForPaper\HALO AI output\Tumor classification\';

files = dir([path '*.rois']);

for i = 1:size(files,1)
    roiList = parseRegionOfInterestList([ files(i).folder filesep files(i).name]);
    for ii = 1:4
        roi = roiList.get(ii);
        mkdir(['X:\PDAC Combo\az neg desi\dpo\rois\' files(i).name(1:end-5) ' ' roi.name])
        cd(['X:\PDAC Combo\az neg desi\dpo\rois\' files(i).name(1:end-5) ' ' roi.name])
        save('roi','roi')
    end
end