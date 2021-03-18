
MyFoldersInfo = dir('*');
MyFoldersInfo = MyFoldersInfo(3:end,:);

for i = 1:size(MyFoldersInfo,1)
    
    cd([MyFoldersInfo(i).folder filesep MyFoldersInfo(i).name])
    
    load('roi')
    
    % 1st quadrant
    
    roi_mask = zeros(roi.width,roi.height);
    roi_mask(1:roi.width/2,1:roi.height/2) = 1;
    
    roi = RegionOfInterest(roi.width,roi.height);
    roi.addPixels(roi_mask)
    
    figure; imagesc(roi.pixelSelection); axis image off
    
    mkdir([MyFoldersInfo(i).folder filesep MyFoldersInfo(i).name '-1'])
    cd([MyFoldersInfo(i).folder filesep MyFoldersInfo(i).name '-1'])
    
    save('roi','roi')
    
    % 2nd quadrant
    
    roi_mask = zeros(roi.width,roi.height);
    roi_mask(1:roi.width/2,roi.height/2+1:end) = 1;
    
    roi = RegionOfInterest(roi.width,roi.height);
    roi.addPixels(roi_mask)
    
    figure; imagesc(roi.pixelSelection); axis image off
    
    mkdir([MyFoldersInfo(i).folder filesep MyFoldersInfo(i).name '-2'])
    cd([MyFoldersInfo(i).folder filesep MyFoldersInfo(i).name '-2'])
    
    save('roi','roi')
    
    % 3rd quadrant
    
    roi_mask = zeros(roi.width,roi.height);
    roi_mask(roi.width/2+1:end,1:roi.height/2) = 1;
    
    roi = RegionOfInterest(roi.width,roi.height);
    roi.addPixels(roi_mask)
    
    figure; imagesc(roi.pixelSelection); axis image off
    
    mkdir([MyFoldersInfo(i).folder filesep MyFoldersInfo(i).name '-3'])
    cd([MyFoldersInfo(i).folder filesep MyFoldersInfo(i).name '-3'])
    
    save('roi','roi')
    
    % 4nd quadrant
    
    roi_mask = zeros(roi.width,roi.height);
    roi_mask(roi.width/2+1:end,roi.height/2+1:end) = 1;
    
    roi = RegionOfInterest(roi.width,roi.height);
    roi.addPixels(roi_mask)
    
    figure; imagesc(roi.pixelSelection); axis image off
    
    mkdir([MyFoldersInfo(i).folder filesep MyFoldersInfo(i).name '-4'])
    cd([MyFoldersInfo(i).folder filesep MyFoldersInfo(i).name '-4'])
    
    save('roi','roi')
    
end