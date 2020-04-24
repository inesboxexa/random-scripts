
D = 'X:\Beatson\SLC7a5 study\AZ data\13C glutamine labelled\120um dpo\single ion images\testrow01NPLconvert\tissue only\no norm\U13C-glutamine SLC7a5';
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
fileID = fopen('quality_scores.txt','w');
fprintf(fileID,txt_row, table');
fclose(fileID);

%%

scores_file = 'X:\Beatson\SLC7a5 study\AZ data\13C glutamine labelled\120um dpo\single ion images\testrow01NPLconvert\tissue only\no norm\U13C-glutamine SLC7a5\quality_scores_final.txt';
scores_struct = tdfread(scores_file,'\t');

mean_intensity_table = 'X:\Beatson\SLC7a5 study\Rafas paper rebuttal\original\mean_intensities_aminoacids_az120umdata_no_norm_apc_apc-kras_unscored_table.xlsx';
[num,txt,raw] = xlsread(mean_intensity_table);

scores_indexes = strcmpi(string(scores_struct.Final_Score),"Good");
scores_meas_mz = scores_struct.meas_mz;

[~, table_row_indexes] = ismembertol(scores_meas_mz(scores_indexes),num(:,1),1e-7);

col_num = size(raw,2)+1;
for rowi = table_row_indexes'
    raw{rowi+1,col_num} = 'Good'; % Correct for the fact the first row are the column names!
end
raw{1,col_num} = 'Quality';

cd('X:\Beatson\SLC7a5 study\Rafas paper rebuttal\original\')
xlswrite('mean_intensities_aminoacids_az120umdata_no_norm_apc_apc-kras_scored_table.xlsx',raw)


