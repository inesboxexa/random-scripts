dataPath = 'D:\veal-brain-homogenate-study\data\';
WatersConvert = '"X:\WatersConvert\watersconvert.exe"';
addpath(genpath('X:\WatersConvert\'));

rawToProcess = dir([dataPath '*.raw']); %find all raw files to convert

for i = 1:size(rawToProcess)
    rawName = [ '"' dataPath rawToProcess(i).name '"' ]; %create loaction of RAW file
    system([WatersConvert ' ' rawName]) %call imzML converter through command line
    disp(['! file iter: ' num2str(i)])
end
    