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

% data folder
d_fld = 'imap_lo';
% CDF  folder
c_fld = 'matlab_cdf391_patch-arm64';

addpath([dir d_fld]);
addpath([dir c_fld]);

%==============================================
%% Import & Process CDFs

% DE
[l1b_dede_data, l1b_dede_info] = ...
    spdfcdfread('imap_lo_l1b_de_20260121-repoint00133_v003.cdf');
l1b_dede_var = l1b_dede_info.Variables(:, [1 4]);

% DE rate
[l1b_dera_data, l1b_dera_info] = ...
    spdfcdfread('imap_lo_l1b_derates_20260121-repoint00133_v001.cdf');
l1b_dera_var = l1b_dera_info.Variables(:, [1 4]);

% histrates
[l1b_hist_data, l1b_hist_info] = ...
    spdfcdfread('imap_lo_l1b_histrates_20260121-repoint00133_v001.cdf');
l1b_hist_var = l1b_hist_info.Variables(:, [1 4]);

% INST state vector
[l1b_inst_data, l1b_inst_info] = ...
    spdfcdfread('imap_lo_l1b_instrument-status-summary_20260121-repoint00133_v001.cdf');
l1b_inst_var = l1b_inst_info.Variables(:, [1 4]);

% monitor rates
[l1b_moni_data, l1b_moni_info] = ...
    spdfcdfread('imap_lo_l1b_monitorrates_20260121-repoint00133_v001.cdf');
l1b_moni_var = l1b_moni_info.Variables(:, [1 4]);

% NHK
[l1b_inhk_data, l1b_inhk_info] = ...
    spdfcdfread('imap_lo_l1b_nhk_20260121-repoint00133_v001.cdf');
l1b_inhk_var = l1b_inhk_info.Variables(:, [1 4]);

% SHK
[l1b_ishk_data, l1b_ishk_info] = ...
    spdfcdfread('imap_lo_l1b_shk_20260121-repoint00133_v001.cdf');
l1b_ishk_var = l1b_ishk_info.Variables(:, [1 4]);

%%
%>>> INCOMPLETE BELOW, NOT IN USE

l1a_hist_shcoarse   = cell2mat(l1a_hist_data( 8));
l1a_hist_count      = cell2mat(l1a_hist_data( 9));
l1a_hist_chksum     = cell2mat(l1a_hist_data(10));
l1a_hist_datadata   = cell2mat(l1a_hist_data(11));
l1a_hist_epoch      = cell2mat(l1a_hist_data(12));
l1a_hist_data_index = cell2mat(l1a_hist_data(13));

l1b_hist_spin_ang_bin = cell2mat(l1b_hist_data(1));
l1b_hist_met          = cell2mat(l1b_hist_data(2));
l1b_hist_data_ave_amp = cell2mat(l1b_hist_data(3));
l1b_hist_data_ave_amp(abs(l1b_hist_data_ave_amp) > 1e30) = NaN; % fill values
l1b_hist_cunt_per_bin = cell2mat(l1b_hist_data(4));
l1b_hist_epoch        = cell2mat(l1b_hist_data(5));
l1b_hist_data_spn_ang = cell2mat(l1b_hist_data(6));


%==============================================
