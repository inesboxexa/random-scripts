
D = 'X:\Beatson\SLC7a5 study\AZ data\13C glutamine labelled\120um dpo\single ion images\testrow01NPLconvert\tissue only\RMS\U13C-glutamine SLC7a5';
S = dir(fullfile(D,'*.png')); % pattern to match filenames.

table = string([]);
for k = 1:numel(S)
    F = fullfile(D,S(k).name);
    I = imread(F);
    imshow(I)
    quality_score = input('Quality: ','s');
    underscore_indices = strfind(S(k).name,'_');
    table = [ 
        table 
        num2str(S(k).name(1:underscore_indices(1)-1),'%1.12f') string(quality_score)
        ];
end
cd(D)
txt_row = strcat(repmat('%s\t',1,size(table,2)-1),'%s\n');
fileID = fopen('quality_scores_2.txt','w');
fprintf(fileID,txt_row, table');
fclose(fileID);

%%

% Scores

scores_file = 'X:\Beatson\SLC7a5 study\AZ data\13C glutamine labelled\120um dpo\single ion images\testrow01NPLconvert\tissue only\RMS\U13C-glutamine SLC7a5\quality_scores_final.txt';
scores_struct = tdfread(scores_file,'\t');
scores_indexes = strcmpi(string(scores_struct.Final_Score),"Good");
scores_meas_mz = scores_struct.meas_mz;

% Mean Intensity Table

mean_intensity_table = 'X:\Beatson\SLC7a5 study\AZ data\13C glutamine labelled\120um dpo\ttest\tissue only\RMS\ttest small intestine vs colon draft.xlsx';
[num,txt,raw] = xlsread(mean_intensity_table);
% mean_intensity_struct = tdfread(mean_intensity_table,'\t');

% [~, table_row_indexes] = ismembertol(scores_meas_mz(scores_indexes),num(:,1),1e-10);
% if size(unique(scores_meas_mz(scores_indexes)),1)~=size(unique(num(table_row_indexes,1)),1); disp('!!! m/z values missing'); end

col_num = size(raw,2)+1;
for rowi = 2:size(raw,1)
    if ismembertol(num(rowi-1,1),scores_meas_mz(scores_indexes),1e-10) % Meas m/z Comparison
        raw{rowi,col_num} = 'Good'; % Correct for the fact the first row are the column names!
    end
end
raw{1,col_num} = 'Quality';

cd('X:\Beatson\SLC7a5 study\Rafas paper rebuttal\original\')
xlswrite('ttest small intestine vs colon draft scored.xlsx',raw)
