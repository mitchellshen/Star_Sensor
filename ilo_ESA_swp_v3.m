% ============================
% == ver. 2024-01-10 by MMS ==
% ============================

%% Initial setup
clc; clear all; close all; format long;

% IMPORT
dirloc = '/Users/mitchellshen/Downloads/';

for k = 2:3
    switch k
        case 1
    inst = readtable([ dirloc 'Instrument_Hybrid_T102_R034_ILO_APP_NHK_20240406T211731_.csv']);
    snif = readtable([ dirloc 'Sniffer_Hybrid_T102_R034_ILO_RAW_CNT_20240406T211732_.csv']);
        case 2
    inst = readtable([ dirloc 'Instrument_Hybrid_T102_R068_ILO_APP_NHK_20240413T221448_.csv']);
    snif = readtable([ dirloc 'Sniffer_Hybrid_T102_R068_ILO_RAW_CNT_20240413T221437_.csv']);
        case 3
    inst = readtable([ dirloc 'Instrument_Hybrid_T102_R071_ILO_APP_NHK_20240414T160433_.csv']);
    snif = readtable([ dirloc 'Sniffer_Hybrid_T102_R071_ILO_RAW_CNT_20240414T160422_.csv']);    
        case 4
    inst = readtable([ dirloc 'Instrument_Hybrid_T102_R072_ILO_APP_NHK_20240414T171310_.csv']);
    snif = readtable([ dirloc 'Sniffer_Hybrid_T102_R072_ILO_RAW_CNT_20240414T171259_.csv']);            
        case 5
    inst = readtable([ dirloc 'Instrument_Hybrid_T102_R073_ILO_APP_NHK_20240414T182226_.csv']);
    snif = readtable([ dirloc 'Sniffer_Hybrid_T102_R073_ILO_RAW_CNT_20240414T182216_.csv']);            
    end

SHCOARSE_INST   = table2array(inst(:,  1));
BHV_ESA_NEG_DAC = table2array(inst(:,111));
BHV_ESA_POS_DAC = table2array(inst(:,112));

SHCOARSE_SNIF   = table2array(snif(:,  1));
SILVER_TRIPLE_0 = table2array(snif(:, 27));
SILVER_TRIPLE_1 = table2array(snif(:, 28));
SILVER_TRIPLE_2 = table2array(snif(:, 29));
SILVER_TRIPLE_3 = table2array(snif(:, 30));
SILVER_TRIPLE   = SILVER_TRIPLE_0+SILVER_TRIPLE_1+SILVER_TRIPLE_2+SILVER_TRIPLE_3;

% PROCESS
[idx_U_NEG,U_NEG_range] = findgroups(BHV_ESA_NEG_DAC);
[idx_U_POS,U_POS_range] = findgroups(BHV_ESA_POS_DAC);

idx_U_NEG_num = 1:1:10;
    switch k
        case 1
    idx_U_POS_num = [1 2 4 5  7  8 10 11 13 14 16 17 19 20 22:25];
        case 2
    idx_U_POS_num = [1 3 6 8 11 12 14 15 17 18 20 21 23 24 26:29];
        case 3
    idx_U_POS_num = [1 3 6 8 11 13 16 18 21 23 26 28 31 33 36 38 40 41];
        case 4
    idx_U_POS_num = [1 3 6 8 11 13 16 18 21 23 26 28 31 33 36 38 40 41];
        case 5
    idx_U_POS_num = [1 3 6 8 11 13 16 18 21 23 26 28 31 33 36 38 40 41];    
    end
    
U_NEG_range = U_NEG_range(idx_U_NEG_num);
U_POS_range = U_POS_range(idx_U_POS_num);
idx_time_l  = zeros(length(idx_U_NEG_num),length(idx_U_POS_num));
idx_time_r  = zeros(length(idx_U_NEG_num),length(idx_U_POS_num));
bin_cnts_Q0 = zeros(length(idx_U_NEG_num),length(idx_U_POS_num));
bin_cnts_all= zeros(length(idx_U_NEG_num),length(idx_U_POS_num));

for j = 1:length(idx_U_POS_num)
for i = 1:length(idx_U_NEG_num) 
    tmp = find(idx_U_NEG==idx_U_NEG_num(i) & idx_U_POS==idx_U_POS_num(j));
    idx_time_l(i,j) = SHCOARSE_INST(tmp(1));
    idx_time_r(i,j) = SHCOARSE_INST(tmp(end));
    U_NEG(i,j) = mean(BHV_ESA_NEG_DAC(tmp(1):tmp(end)));
    U_POS(i,j) = mean(BHV_ESA_POS_DAC(tmp(1):tmp(end)));
    clear tmp
    tmp = find(SHCOARSE_SNIF >= idx_time_l(i,j) & SHCOARSE_SNIF <= idx_time_r(i,j));
    bin_cnts_Q0(i,j) = mean(SILVER_TRIPLE_0(tmp(1):tmp(end)));
    bin_cnts_all(i,j) = mean(SILVER_TRIPLE(tmp(1):tmp(end)));
    clear tmp
end
end

bin_cnts_proc = bin_cnts_all;

% ESA voltage
U_NEG_matrix = ones(1,size(bin_cnts_proc,2)) .* U_NEG_range;
U_NEG_collapse = reshape(U_NEG_matrix,1,[]);

% outer rotational angle
U_POS_matrix = ones(1,size(bin_cnts_proc,1)) .* U_POS_range;
U_POS_collapse = reshape(U_POS_matrix',1,[]);
I_cnt_collapse = reshape(bin_cnts_proc,1,[]);


% PLOT
% close all
dotsize     = 800;   % For passband purposes. nominal = 60 & extended = 20

figure(1); set(gcf,'Position', [50, 150, 1100, 500 ]);

ax(k) = subplot(1,2,k);
scatter(U_NEG_collapse,U_POS_collapse,dotsize,I_cnt_collapse,...
        'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'marker','s')
set(gca,'color',0*[1 1 1]); hold on; 
box on; axis equal; %alpha(.9)
set(gca,'lineWidth',1.0,'FontSize',13,'FontWeight','bold');
hcb = colorbar; 
caxis([0 ceil(max(max(bin_cnts_proc))/10)*10]);
colormap(ax(k),turbo); 
ylabel(hcb,'Binned Summed Silver Triple Counts','FontSize',13)
xlabel('U- (V)');
ylabel('U+ (V)');
xlim([min(U_NEG_range)-100 max(U_NEG_range)+100]); 
ylim([min(U_POS_range)-100 max(U_POS_range)+100]); 
title('U+/U- Sweep (930H Beam)')
set(gca, 'XDir','normal'); set(gca, 'YDir','normal'); %reverse
a = gca; b = copyobj(a, gcf);
set(b,'Xcolor',1*[1 1 1],'YColor',1*[1 1 1],'XTickLabel',[],'YTickLabel',[],...
          'XLabel',[],'YLabel',[],'Title',[],'TickLength',[0.02 0.05]);

switch k
    case 1
        before = I_cnt_collapse;
        before_norm = before./max(max(before));
    case 2
        after  = I_cnt_collapse;
        after_norm = after./max(max(after));
end

% ax(2) = subplot(1,2,2);
% pcolor(U_NEG(:,end)-50,U_POS(end,:)-50,bin_cnts_all')
% shading flat;
% box on; axis equal;
% set(gca,'color',0*[1 1 1]); hold on; 
% xlabel('U- (V)')
% ylabel('U+ (V)')
% xlim([min(U_NEG_range)-100 max(U_NEG_range)+100]); 
% ylim([min(U_POS_range)-100 max(U_POS_range)+100]); 
% set(gca,'lineWidth',1.0,'FontSize',13,'FontWeight','bold');
% title('U+/U- Sweep (930H Beam)')
% % title('930H Beam, U+/U- Sweep from eStep 6-7 Sum4Quads')
% hcb = colorbar; caxis([0 ceil(max(max(bin_cnts_all))/10)*10]);
% colormap(ax(2),turbo); 
% ylabel(hcb,'Binned Summed Silver Triple Counts','FontSize',13)
% a = gca; b = copyobj(a, gcf);
% set(b,'Xcolor',1*[1 1 1],'YColor',1*[1 1 1],'XTickLabel',[],'YTickLabel',[],...
%           'XLabel',[],'YLabel',[],'Title',[],'TickLength',[0.02 0.05]);

end
% 
% ax(3) = subplot(1,3,3);
% scatter(U_NEG_collapse,U_POS_collapse,dotsize,after_norm./before_norm-1,...
%         'filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'marker','s')
% set(gca,'color',0*[1 1 1]); hold on; 
% box on; axis equal; %alpha(.9)
% set(gca,'lineWidth',1.0,'FontSize',13,'FontWeight','bold');
% hcb = colorbar; 
% % caxis([0 ceil(max(max(bin_cnts_all))/10)*10]);
% colormap(ax(3),bluewhitered); 
% ylabel(hcb,'Binned Summed Silver Triple Counts','FontSize',13)
% xlabel('U- (V)');
% ylabel('U+ (V)');
% xlim([min(U_NEG_range)-100 max(U_NEG_range)+100]); 
% ylim([min(U_POS_range)-100 max(U_POS_range)+100]); 
% title('U+/U- Sweep (930H Beam)')
% set(gca, 'XDir','normal'); set(gca, 'YDir','normal'); %reverse
% set(gca,'ColorScale','log')
% a = gca; b = copyobj(a, gcf);
% set(b,'Xcolor',1*[1 1 1],'YColor',1*[1 1 1],'XTickLabel',[],'YTickLabel',[],...
%           'XLabel',[],'YLabel',[],'Title',[],'TickLength',[0.02 0.05]);
