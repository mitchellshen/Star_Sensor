close all; clearvars; clc; format long; warning('off','all')

%% Set directory
computer = 'M';

if computer == 'M'      % personal MAC
    %dir = '/Users/mitchellshen/Dropbox (Personal)/Dust_experiment/';
    dir = '/Users/mitchellshen/Library/CloudStorage/Dropbox-Personal/Dust_experiment/';
elseif computer == 'W'  % personal WIN
    dir = 'C:\Users\miche\Dropbox (個人)\Dust_experiment\';
elseif computer == 'L'  % Lab WIN
    dir = 'C:\Users\ms3648\Dropbox (Personal)\Dust_experiment\';
end

d_fld = 'imap_lo';
c_fld = 'matlab_cdf391_patch-arm64';

addpath([dir d_fld]);
addpath([dir c_fld]);

%==============================================
%% Import CDFs
[l1b_inst_data, l1b_inst_info] = ...
    spdfcdfread('imap_lo_l1b_instrument-status-summary_20260121-repoint00133_v001.cdf');
    % spdfcdfread('imap_lo_l1b_instrument-status-summary_20251104-repoint00038_v001.cdf');

[l1a_nhk_data, l1a_nhk_info] = ...
    spdfcdfread('imap_lo_l1a_nhk_20260121-repoint00133_v001.cdf');
    % spdfcdfread('imap_lo_l1a_nhk_20251104-repoint00038_v001.cdf');

[l1a_shk_data, l1a_shk_info] = ...
    spdfcdfread('imap_lo_l1a_shk_20260121-repoint00133_v001.cdf');
    % spdfcdfread('imap_lo_l1a_shk_20251104-repoint00038_v001.cdf');

[l1b_nhk_data, l1b_nhk_info] = ...
    spdfcdfread('imap_lo_l1b_nhk_20260121-repoint00133_v001.cdf');
    % spdfcdfread('imap_lo_l1b_nhk_20251104-repoint00038_v001.cdf');

[l1b_shk_data, l1b_shk_info] = ...
    spdfcdfread('imap_lo_l1b_shk_20260121-repoint00133_v001.cdf');
    % spdfcdfread('imap_lo_l1b_shk_20251104-repoint00038_v001.cdf');


l1b_inst_var = l1b_inst_info.Variables;
% l1b_inst_info.Variables(:, [1 4]);

l1a_nhk_var = l1a_nhk_info.Variables(:, [1 4]);
l1a_shk_var = l1a_shk_info.Variables(:, [1 4]);
l1b_nhk_var = l1b_nhk_info.Variables(:, [1 4]);
l1b_shk_var = l1b_shk_info.Variables(:, [1 4]);

%==============================================
%% TOF threshold
figure(1)
% -- Initial check --
% varnum = [1 2 6 8];
% varnum_l1a = [273 271 270 272];
% varname = {'tof0_thr', 'tof2_thr', 'tof3_thr', 'tof1_thr' };

% -- Flight check -- 
varnum = [19 6 2 14];
varnum_l1a = [280 279 278 277];
varname = {'tof0_thr', 'tof1_thr', 'tof2_thr', 'tof3_thr' };

for i = 1:length(varnum)
subplot(4,1,i)
plot( cell2mat(l1b_inst_data(varnum(i)))    ,'LineWidth',2); hold on;
plot( cell2mat(l1a_nhk_data(varnum_l1a(i))) ,'LineWidth',2); hold on;
ylabel(varname(i), 'Interpreter', 'none'); grid on;
legend({'l1b (after)','l1a (before)'});
end

% Validated with issues found on 2026-01-26
% TOF threshold (0,1,2,3) = (15,12,7,11)

%==============================================
%% Anode threshold
figure(2)
% -- Initial check --
% varnum = [5 13 15 19];
% varnum_l1a = [251 249 250 252];
% varname = {'an_b3_thr', 'an_a_thr', 'an_b0_thr', 'an_c_thr' };

% -- Flight check --
varnum = [16 17 8 15];
varnum_l1a = [256 259 257 258];
varname = {'an_a_thr', 'an_c_thr', 'an_b0_thr', 'an_b3_thr' };


for i = 1:length(varnum)
subplot(4,1,i)
plot( cell2mat(l1b_inst_data(varnum(i)))    ,'LineWidth',1); hold on;
plot( cell2mat(l1a_nhk_data(varnum_l1a(i))) ,'LineWidth',2, 'linestyle','--'); hold on;
ylabel(varname(i), 'Interpreter', 'none'); grid on;
legend({'l1b (after)','l1a (before)'});
end

% Validated on 2026-01-26
% Anode threshold (0,1,2,3) = (0,5,0,0)

%==============================================
%% MCP PAC voltages
figure(3)
% -- Initial check --
% varnum = [4 9 14];
% varnum_l1a = [197 195 223 242];
% varname = {'mcp_vset', 'pac_vset', 'mcp_v', 'tof_mcp_v' };

% -- Flight check --
varnum = [20 7 10];
varnum_l1a = [204 202 230 249];
varname = {'mcp_vset', 'pac_vset', 'mcp_v', 'tof_mcp_v' };


for i = 1:length(varnum)
subplot(3,1,i)
plot( cell2mat(l1b_inst_data(varnum(i)))    ,'LineWidth',2); hold on;
plot( cell2mat(l1a_nhk_data(varnum_l1a(i))) ,'LineWidth',2, 'linestyle','--'); hold on;
ylabel(varname(i), 'Interpreter', 'none'); grid on;
legend({'l1b (after)','l1a (before)'});
end
subplot(3,1,3)
plot( cell2mat(l1a_nhk_data(varnum_l1a(4))) ,'LineWidth',2, 'linestyle','--'); hold on;
ylabel(varname(4), 'Interpreter', 'none'); grid on;
legend({'l1b (after)','l1a (mcp_v)','l1a (tof_mcp_v)'}, 'Interpreter', 'none');

% Validated on 2026-01-26
% MCP VSET 3144, PAC 12kV, TOF_MCP_V ~2300

%==============================================
%% DEF and PMT
figure(4)
% -- Initial check --
% varnum = [10 11 18];
% varnum_l1a = [96 98 94];
% varname = {'bhv_def_pos_dac', 'bhv_pmt_dac', 'bhv_def_neg_dac'};

% -- Flight check --
varnum = [13 18 5];
varnum_l1a = [103 101 105];
varname = {'bhv_def_pos_dac', 'bhv_def_neg_dac', 'bhv_pmt_dac'};

for i = 1:length(varnum)
subplot(3,1,i)
plot( cell2mat(l1b_inst_data(varnum(i)))    ,'LineWidth',2); hold on;
plot( cell2mat(l1a_nhk_data(varnum_l1a(i))) ,'LineWidth',2); hold on;
ylabel(varname(i), 'Interpreter', 'none'); grid on;
legend({'l1b (after)','l1a (before)'});
end


%==============================================
%% IFB & PCC
figure(5)
% -- Initial check --
% varnum = [17 21 22];
% varnum_l1a = [203 280 282];
% varname = {'ifb_hot_spot_t', 'coarse_pot_pri', 'fine_pot_pri'};

% -- Flight check --
varnum = [3 21 22];
varnum_l1a = [210 287 289];
varname = {'ifb_hot_spot_t', 'coarse_pot_pri', 'fine_pot_pri'};

for i = 1:length(varnum)
subplot(3,1,i)
plot( cell2mat(l1b_inst_data(varnum(i)))    ,'LineWidth',2); hold on;
plot( cell2mat(l1a_nhk_data(varnum_l1a(i))) ,'LineWidth',2); hold on;
ylabel(varname(i), 'Interpreter', 'none'); grid on;
legend({'l1b (after)','l1a (before)'});
end

%==============================================
%% Char test
% -- Initial check --
% eng_lut = cell2mat(l1b_inst_data( 3));
% sci_lut = cell2mat(l1b_inst_data(12));
% op_mode = cell2mat(l1b_inst_data(16));
% fsw_ver = cell2mat(l1b_inst_data(20));

% -- Flight check --
eng_lut = cell2mat(l1b_inst_data( 9));
sci_lut = cell2mat(l1b_inst_data(12));
op_mode = cell2mat(l1b_inst_data( 1));
fsw_ver = cell2mat(l1b_inst_data(11));

%==============================================
%% shcoarse and epoch
figure(6)
% -- Initial check --
% varnum = [7 23];
% varnum_l1a = [1 354];
% varname = {'shcoarse', 'epoch'};

% -- Flight check --
varnum = [4 23];
varnum_l1a = [8 354];
varname = {'shcoarse', 'epoch'};

for i = 1:length(varnum)
subplot(2,1,i)
plot( cell2mat(l1b_inst_data(varnum(i)))    ,'LineWidth',2); hold on;
plot( cell2mat(l1a_nhk_data(varnum_l1a(i))) ,'LineWidth',2); hold on;
ylabel(varname(i), 'Interpreter', 'none'); grid on;
legend({'l1b (after)','l1a (before)'});
end

%% New validation

[l1b_inst_data_N, l1b_inst_info_N] = ...
    spdfcdfread('imap_lo_l1b_instrument-status-summary_20251110-repoint00044_v001.cdf');
