
% Recovering the amplitude of the peaks saved in the datacube

cd('X:\ICR Breast PDX\Negative Analysis\Data Processing Outputs CA Pellets and PDX\spectra details\')

files_struct = dir;

% pdx-1282 / pdx-1458 / cell-pellet-1282 / cell-pellet-pdx-1458

new_file_indexs = [ 4 8 16 12 6 10 18 14 3 7 15 11 5 9 17 13 ];
y = [];
read_x = 1;
px = [];
py = [];
aux_index = 0;

for file_index = new_file_indexs
    
    aux_index = aux_index + 1;
       
    if read_x
        
        load([ files_struct(file_index).folder filesep files_struct(file_index).name '/tissue only/totalSpectrum_mzvalues']);
        x = totalSpectrum_mzvalues';
        
    end
    
    load([ files_struct(file_index).folder filesep files_struct(file_index).name '/tissue only/totalSpectrum_intensities']);
    load([ files_struct(file_index).folder filesep files_struct(file_index).name '/tissue only/pixels_num']);
    
    y0 = totalSpectrum_intensities'./pixels_num;
    
    y = [ y, y0 ];
    
    load([ files_struct(file_index).folder filesep files_struct(file_index).name '/tissue only/datacubeonly_peakDetails']);
    
    for peak_index = 1:size(datacubeonly_peakDetails,1)
              
        pmz = datacubeonly_peakDetails(peak_index,2);
        
        [~, peak_index0] = min(abs(x-pmz));
        
        if read_x
            
            px = datacubeonly_peakDetails(:,2);
            read_x = 0;
            
        end
        
        py(peak_index,aux_index) = y0(peak_index0);
        
    end
    
end

y01 = (y-min(y,[],1))./max(y-min(y,[],1),[],1);

%%

figure;
subplot(2,1,1)
plot(x,[y(:,1:4),-y(:,5:8)]); hold on; plot(px,py(:,1:4),'x');
subplot(2,1,2)
plot(x,[y01(:,1:4),-y01(:,5:8)])

%% Anova (2-way)

py_pdx_only = py(:,1:8);

g1 = { '1282'; '1282'; '1282'; '1282'; '1458'; '1458'; '1458'; '1458' }; % genetic model
g2 = { 's1'; 's2'; 's3'; 's4'; 's1'; 's2'; 's3'; 's4' }; % session

for peak_index = 1:size(py_pdx_only,1)

    [ p0, tbl0, stats0, terms0 ] = anovan(py_pdx_only(peak_index,:)', {g1, g2}, 'display','off');

    anova_pvalue(peak_index,1:2) = p0;

end

%%

anova_index_1 = logical((anova_pvalue(:,1)<0.05).*(anova_pvalue(:,2)>0.05)); % peaks that distinguish genetic model
anova_index_2 = logical((anova_pvalue(:,2)<0.05).*(anova_pvalue(:,1)>0.05)); % peaks that distinguish sessions

figure;
subplot(2,1,1)
plot(x,[y(:,1:4),-y(:,5:8)]); hold on; plot(px(anova_index_1),py(anova_index_1,1:4),'x'); plot(px(anova_index_2),py(anova_index_2,1:4),'o');
subplot(2,1,2)
plot(x,[y01(:,1:4),-y01(:,5:8)])

peaks2keep = datacubeonly_peakDetails(~anova_index_2,2);




