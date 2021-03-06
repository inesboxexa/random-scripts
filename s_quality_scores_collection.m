
D = 'X:\Beatson\SLC7a5 study\AZ data\APC_KRAS vs APC_KRAS_SLC7a5\DPO\single ion images\AZ APC-KRAS-SLC7a5 Colon\tissue only\RMS\Amino Acids';
S = dir(fullfile(D,'*.png')); % pattern to match filenames.

table = ["meas_mz","Final_Score"];
for k = 1:numel(S)     
    underscore_indices = strfind(S(k).name,'_');
    if ~isnan(str2double(S(k).name(1:underscore_indices(1)-1)))
        F = fullfile(D,S(k).name);
        I = imread(F);
        imshow(I)
        quality_score = input('Quality: ','s');
        table = [
            table
            num2str(S(k).name(1:underscore_indices(1)-1),'%1.12f') string(quality_score)
            ];
    end
end
cd(D)
txt_row = strcat(repmat('%s\t',1,size(table,2)-1),'%s\n');
fileID = fopen('quality_scores_2.txt','w');
fprintf(fileID,txt_row, table');
fclose(fileID);

%%

% Scores

scores_file = 'X:\Beatson\SLC7a5 study\AZ data\APC_KRAS vs APC_KRAS_SLC7a5\DPO\single ion images\AZ APC-KRAS-SLC7a5 Colon\tissue only\RMS\Amino Acids\quality_scores_final.txt';
scores_struct = tdfread(scores_file,'\t');
scores_indexes = strcmpi(string(scores_struct.Final_Score),"G");
scores_meas_mz = scores_struct.meas_mz(scores_indexes,:);

% scores_meas_mz = [];
% for scores_i = 1:length(scores_meas_mz0)
%     scores_meas_mz(scores_i,1) = str2double(scores_meas_mz0(scores_i,:));
% end

% Mean Intensity Table

mean_intensity_table = 'X:\Beatson\SLC7a5 study\AZ data\APC_KRAS vs APC_KRAS_SLC7a5\DPO\ttest\tissue only\RMS\ttest APC-KRAS-SLC7a5 Colon vs APC-KRAS Colon.xlsx';
[num,txt,raw] = xlsread(mean_intensity_table);
% mean_intensity_struct = tdfread(mean_intensity_table,'\t');

% [~, table_row_indexes] = ismembertol(scores_meas_mz(scores_indexes),num(:,1),1e-10);
% if size(unique(scores_meas_mz(scores_indexes)),1)~=size(unique(num(table_row_indexes,1)),1); disp('!!! m/z values missing'); end
i = 0;
col_num = size(raw,2)+1;
for rowi = 2:size(raw,1)
    if ismembertol(num(rowi-1,1),scores_meas_mz,1e-6) % Meas m/z Comparison
        raw{rowi,col_num} = 'Good'; % Correct for the fact the first row are the column names!
        i = i+1; disp(i); 
    end
end
raw{1,col_num} = 'Quality';

cd('X:\Beatson\SLC7a5 study\AZ data\APC_KRAS vs APC_KRAS_SLC7a5\DPO\ttest\tissue only\RMS\')
xlswrite('ttest APC-KRAS-SLC7a5 Colon vs APC-KRAS Colon scored.xlsx',raw)
