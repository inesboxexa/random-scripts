
load('X:\ICR Breast PDX\Data Processing Outputs\all pdx pilot\uva\tissue only\no norm\study 4 - neg desi - wt vs mutant - all combi.mat')
nonorm_data = double(table(2:end,1:10));

load('X:\ICR Breast PDX\Data Processing Outputs\all pdx pilot\uva\tissue only\tic\study 4 - neg desi - wt vs mutant - all combi.mat')
tic_data = double(table(2:end,1:10));

load('X:\ICR Breast PDX\Data Processing Outputs\all pdx pilot\uva\tissue only\RMS\study 4 - neg desi - wt vs mutant - all combi.mat')
rms_data = double(table(2:end,1:10));

load('X:\ICR Breast PDX\Data Processing Outputs\all pdx pilot\uva\tissue only\pqn median\study 4 - neg desi - wt vs mutant - all combi.mat')
pqnmedian_data = double(table(2:end,1:10));

load('X:\ICR Breast PDX\Data Processing Outputs\all pdx pilot\uva\tissue only\pqn mean\study 4 - neg desi - wt vs mutant - all combi.mat')
pqnmean_data = double(table(2:end,1:10));

load('X:\ICR Breast PDX\Data Processing Outputs\all pdx pilot\uva\tissue only\z-score\study 4 - neg desi - wt vs mutant - all combi.mat')
z_data = double(table(2:end,1:10));

load('X:\ICR Breast PDX\Data Processing Outputs\all pdx pilot\uva\tissue only\log\study 4 - neg desi - wt vs mutant - all combi.mat')
log_data = double(table(2:end,1:10));

mz = double(table(2:end,11));

colours = tab10;

figure;
for col = 1:size(nonorm_data,2)
    hold on; 
    plot(mz,nonorm_data(:,col),'.','color',colours(col,:));
end

%%

figure; imagesc(nonorm_data'); colormap(viridis(10))

figure; imagesc(tic_data'); colormap(viridis(10))

%%

xi = 9;
yi = 8;

figure

subplot(2,3,1)
hold on;
plot([0 1],[0 1],'color',[0.75 0.75 0.75]); grid on;
plot(nonorm_data(:,xi),nonorm_data(:,yi),'.k'); grid on;
stem([0.3 0.7],[1.1 1.1],'color',[1 0 0.75]); grid on;
plot([0 1],[0.3 0.3],'color',[1 0 0.75]); grid on;
plot([0 1],[0.7 0.7],'color',[1 0 0.75]); grid on;
ylabel('TR')
xlabel('BR')
axis([0 1 0 1])
title('No normalisation')

subplot(2,3,2) 
hold on;
plot([0 1],[0 1],'color',[0.75 0.75 0.75]); grid on;
plot(tic_data(:,xi),tic_data(:,yi),'.k'); grid on;
stem([0.3 0.7],[1.1 1.1],'color',[1 0 0.75]); grid on;
plot([0 1],[0.3 0.3],'color',[1 0 0.75]); grid on;
plot([0 1],[0.7 0.7],'color',[1 0 0.75]); grid on;
ylabel('TR')
xlabel('BR')
axis([0 1 0 1])
title('TIC')

subplot(2,3,3) 
hold on;
plot([0 1],[0 1],'color',[0.75 0.75 0.75]); grid on;
plot(rms_data(:,xi),rms_data(:,yi),'.k'); grid on;
stem([0.3 0.7],[1.1 1.1],'color',[1 0 0.75]); grid on;
plot([0 1],[0.3 0.3],'color',[1 0 0.75]); grid on;
plot([0 1],[0.7 0.7],'color',[1 0 0.75]); grid on;
ylabel('TR')
xlabel('BR')
axis([0 1 0 1])
title('RMS')

subplot(2,3,4) 
hold on;
plot([0 1],[0 1],'color',[0.75 0.75 0.75]); grid on;
plot(pqnmedian_data(:,xi),pqnmedian_data(:,yi),'.k'); grid on;
stem([0.3 0.7],[1.1 1.1],'color',[1 0 0.75]); grid on;
plot([0 1],[0.3 0.3],'color',[1 0 0.75]); grid on;
plot([0 1],[0.7 0.7],'color',[1 0 0.75]); grid on;
ylabel('TR')
xlabel('BR')
axis([0 1 0 1])
title('pqn median')

subplot(2,3,5) 
hold on;
plot([0 1],[0 1],'color',[0.75 0.75 0.75]); grid on;
plot(pqnmean_data(:,xi),pqnmean_data(:,yi),'.k'); grid on;
stem([0.3 0.7],[1.1 1.1],'color',[1 0 0.75]); grid on;
plot([0 1],[0.3 0.3],'color',[1 0 0.75]); grid on;
plot([0 1],[0.7 0.7],'color',[1 0 0.75]); grid on;
ylabel('TR')
xlabel('BR')
axis([0 1 0 1])
title('pqn mean')

%%

figure

subplot(2,3,1)
hold on;
plot([0 1],[0 1],'color',[0.75 0.75 0.75]); grid on;
plot(nonorm_data(:,1:2),nonorm_data(:,3:4),'.k'); grid on;
stem([0.3 0.7],[1.1 1.1],'color',[1 0 0.75]); grid on;
plot([0 1],[0.3 0.3],'color',[1 0 0.75]); grid on;
plot([0 1],[0.7 0.7],'color',[1 0 0.75]); grid on;
ylabel('TR 1 & 2')
xlabel('TR 3 & 4')
axis([0 1 0 1])
title('No normalisation')

subplot(2,3,2) 
hold on;
plot([0 1],[0 1],'color',[0.75 0.75 0.75]); grid on;
plot(tic_data(:,1:2),tic_data(:,3:4),'.k'); grid on;
stem([0.3 0.7],[1.1 1.1],'color',[1 0 0.75]); grid on;
plot([0 1],[0.3 0.3],'color',[1 0 0.75]); grid on;
plot([0 1],[0.7 0.7],'color',[1 0 0.75]); grid on;
ylabel('TR 1 & 2')
xlabel('TR 3 & 4')
axis([0 1 0 1])
title('TIC')

subplot(2,3,3) 
hold on;
plot([0 1],[0 1],'color',[0.75 0.75 0.75]); grid on;
plot(rms_data(:,1:2),rms_data(:,3:4),'.k'); grid on;
stem([0.3 0.7],[1.1 1.1],'color',[1 0 0.75]); grid on;
plot([0 1],[0.3 0.3],'color',[1 0 0.75]); grid on;
plot([0 1],[0.7 0.7],'color',[1 0 0.75]); grid on;
ylabel('TR 1 & 2')
xlabel('TR 3 & 4')
axis([0 1 0 1])
title('RMS')

subplot(2,3,4) 
hold on;
plot([0 1],[0 1],'color',[0.75 0.75 0.75]); grid on;
plot(pqnmedian_data(:,1:2),pqnmedian_data(:,3:4),'.k'); grid on;
stem([0.3 0.7],[1.1 1.1],'color',[1 0 0.75]); grid on;
plot([0 1],[0.3 0.3],'color',[1 0 0.75]); grid on;
plot([0 1],[0.7 0.7],'color',[1 0 0.75]); grid on;
ylabel('TR 1 & 2')
xlabel('TR 3 & 4')
axis([0 1 0 1])
title('pqn median')

subplot(2,3,5) 
hold on;
plot([0 1],[0 1],'color',[0.75 0.75 0.75]); grid on;
plot(pqnmean_data(:,1:2),pqnmean_data(:,3:4),'.k'); grid on;
stem([0.3 0.7],[1.1 1.1],'color',[1 0 0.75]); grid on;
plot([0 1],[0.3 0.3],'color',[1 0 0.75]); grid on;
plot([0 1],[0.7 0.7],'color',[1 0 0.75]); grid on;
ylabel('TR 1 & 2')
xlabel('TR 3 & 4')
axis([0 1 0 1])
title('pqn mean')

%%

figure

subplot(2,3,1)
hold on;
plot([0 1],[0 1],'color',[0.75 0.75 0.75]); grid on;
plot(nonorm_data(:,5:6),nonorm_data(:,7),'.k'); grid on;
stem([0.3 0.7],[1.1 1.1],'color',[1 0 0.75]); grid on;
plot([0 1],[0.3 0.3],'color',[1 0 0.75]); grid on;
plot([0 1],[0.7 0.7],'color',[1 0 0.75]); grid on;
ylabel('BR 1 & 2')
xlabel('BR 3')
axis([0 1 0 1])
title('No normalisation')

subplot(2,3,2) 
hold on;
plot([0 1],[0 1],'color',[0.75 0.75 0.75]); grid on;
plot(tic_data(:,5:6),tic_data(:,7),'.k'); grid on;
stem([0.3 0.7],[1.1 1.1],'color',[1 0 0.75]); grid on;
plot([0 1],[0.3 0.3],'color',[1 0 0.75]); grid on;
plot([0 1],[0.7 0.7],'color',[1 0 0.75]); grid on;
ylabel('BR 1 & 2')
xlabel('BR 3')
axis([0 1 0 1])
title('TIC')

subplot(2,3,3) 
hold on;
plot([0 1],[0 1],'color',[0.75 0.75 0.75]); grid on;
plot(rms_data(:,5:6),rms_data(:,7),'.k'); grid on;
stem([0.3 0.7],[1.1 1.1],'color',[1 0 0.75]); grid on;
plot([0 1],[0.3 0.3],'color',[1 0 0.75]); grid on;
plot([0 1],[0.7 0.7],'color',[1 0 0.75]); grid on;
ylabel('BR 1 & 2')
xlabel('BR 3')
axis([0 1 0 1])
title('RMS')

subplot(2,3,4) 
hold on;
plot([0 1],[0 1],'color',[0.75 0.75 0.75]); grid on;
plot(pqnmedian_data(:,5:6),pqnmedian_data(:,7),'.k'); grid on;
stem([0.3 0.7],[1.1 1.1],'color',[1 0 0.75]); grid on;
plot([0 1],[0.3 0.3],'color',[1 0 0.75]); grid on;
plot([0 1],[0.7 0.7],'color',[1 0 0.75]); grid on;
ylabel('BR 1 & 2')
xlabel('BR 3')
axis([0 1 0 1])
title('pqn median')

subplot(2,3,5) 
hold on;
plot([0 1],[0 1],'color',[0.75 0.75 0.75]); grid on;
plot(pqnmean_data(:,5:6),pqnmean_data(:,7),'.k'); grid on;
stem([0.3 0.7],[1.1 1.1],'color',[1 0 0.75]); grid on;
plot([0 1],[0.3 0.3],'color',[1 0 0.75]); grid on;
plot([0 1],[0.7 0.7],'color',[1 0 0.75]); grid on;
ylabel('BR 1 & 2')
xlabel('BR 3')
axis([0 1 0 1])
title('pqn mean')


