%% Initial setup
clc; clear all; close all; format long;

% IMPORT
S = dir(('/Volumes/Mitchell/IMAP_Lo_*/**/**/**/*.csv'));
idx = [];
for i = 1:length(S)./10
   tmp = 1+10*(i-1);
   idx = [idx tmp];
end
S = S(idx);
L = {S.folder};
N = {S.name};
X_hidden  = ~cellfun('isempty',strfind(N,'._'));            %#ok<STRCL1>
X_ILO_IFB = ~cellfun('isempty',strfind(N,'ILO_IFB'));       %#ok<STRCL1>
disp(['found ' num2str(sum(X_ILO_IFB),'%04d') ' ILO_IFB files']);
X_ILO_TOF = ~cellfun('isempty',strfind(N,'ILO_TOF'));       %#ok<STRCL1>
disp(['found ' num2str(sum(X_ILO_TOF),'%04d') ' ILO_TOF files']);
X_ILO_NHK = ~cellfun('isempty',strfind(N,'ILO_APP_NHK'));   %#ok<STRCL1>
disp(['found ' num2str(sum(X_ILO_NHK),'%04d') ' ILO_NHK files']);

%% PROCESS
N_ILO_IFB = N(X_ILO_IFB & ~X_hidden); L_ILO_IFB = L(X_ILO_IFB & ~X_hidden);
N_ILO_TOF = N(X_ILO_TOF & ~X_hidden); L_ILO_TOF = L(X_ILO_TOF & ~X_hidden);
N_ILO_NHK = N(X_ILO_NHK & ~X_hidden); L_ILO_NHK = L(X_ILO_NHK & ~X_hidden);

SHCOARSE_IFB = nan(1,1); 
PAC_VSET = nan(1,1); 
MCP_VSET = nan(1,1); 
MCP_VM   = nan(1,1); 
PAC_VM   = nan(1,1); 

for i = 1:108 %length(N_ILO_IFB)
    
ILO = table2array(readtable([char(L_ILO_IFB(end-i+1)) '/' char(N_ILO_IFB(end-i+1))]));
disp(char(N_ILO_IFB(end-i+1)));

SHCOARSE_tmp = ILO(:,1);
PAC_VSET_tmp = ILO(:,38);
MCP_VSET_tmp = ILO(:,40);
MCP_VM_tmp   = ILO(:,65);
PAC_VM_tmp   = ILO(:,68);

SHCOARSE_IFB = [SHCOARSE_tmp;   SHCOARSE_IFB];
PAC_VSET     = [PAC_VSET_tmp;   PAC_VSET];
MCP_VSET     = [MCP_VSET_tmp;   MCP_VSET];
MCP_VM       = [MCP_VM_tmp;     MCP_VM];
PAC_VM       = [PAC_VM_tmp;     PAC_VM];

clear SHCOARSE_tmp PAC_VSET_tmp MCP_VSET_tmp MCP_VM_tmp PAC_VM_tmp
end

UTtime_IFB = datetime(2010,1,1)+SHCOARSE_IFB./(24*3600);


SHCOARSE_TOF = nan(1,1); 
TOF_MCP_CM = nan(1,1); 
TOF_MCP_VM = nan(1,1); 
for i = 1:108 %length(N_ILO_IFB)
    
TOF = table2array(readtable([char(L_ILO_TOF(end-i+1)) '/' char(N_ILO_TOF(end-i+1))]));
disp(char(N_ILO_TOF(end-i+1)));

SHCOARSE_tmp = TOF(:,1);
TOF_MCP_CM_tmp   = TOF(:,24);
TOF_MCP_VM_tmp   = TOF(:,25);

SHCOARSE_TOF = [SHCOARSE_tmp;   SHCOARSE_TOF];
TOF_MCP_CM   = [TOF_MCP_CM_tmp;   TOF_MCP_CM];
TOF_MCP_VM   = [TOF_MCP_VM_tmp;   TOF_MCP_VM];

clear SHCOARSE_tmp TOF_MCP_CM_tmp TOF_MCP_VM_tmp
end

UTtime_TOF = datetime(2010,1,1)+SHCOARSE_TOF./(24*3600);

%%
A = SHCOARSE_IFB;
B = SHCOARSE_TOF;
idx_AB = ismember(A, B, 'rows');
idx_BA = ismember(B, A, 'rows');
AA = A(idx_AB);
BB = B(idx_BA);

UTtime_IFB = UTtime_IFB(idx_BA);
MCP_VSET   = MCP_VSET(idx_BA);
MCP_VM     = MCP_VM(idx_BA);

UTtime_TOF = UTtime_TOF(idx_AB);
TOF_MCP_VM = TOF_MCP_VM(idx_AB);

%%
close all;
figure(1)
plot(UTtime_IFB,MCP_VSET,'^','MarkerSize',10); hold on;
plot(UTtime_IFB,MCP_VM,'+','MarkerSize',10); hold on;
plot(UTtime_TOF,TOF_MCP_VM,'d','MarkerSize',10); hold on;
set(gca,'lineWidth',1,'FontSize',14,'FontWeight','bold');
dl = datetime('22-Jul-2024');
dr = datetime('26-Jul-2024');
interval = days(2);
xtix = dl : interval : dr;
set(gca, 'XTick',xtix)
legend({'MCP_VSET','MCP_VM','TOF_MCP_VM'})
xlim([dl dr])


