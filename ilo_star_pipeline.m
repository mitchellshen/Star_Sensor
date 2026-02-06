close all; clearvars; clc; format long; warning('off','all')

%% Set directory
computer = 'M';

if computer == 'M'      % personal MAC
    %dir = '/Users/mitchellshen/Dropbox (Personal)/Dust_experiment/';
    s_fld = '/Users/mitchellshen/Library/CloudStorage/Dropbox-Personal/Dust_experiment/';
elseif computer == 'W'  % personal WIN
    s_fld = 'C:\Users\miche\Dropbox (個人)\Dust_experiment\';
elseif computer == 'L'  % Lab WIN
    dis_fldr = 'C:\Users\ms3648\Dropbox (Personal)\Dust_experiment\';
end

% data folder
d_fld = ['/Users/mitchellshen/Library/CloudStorage/Dropbox-Personal/'...
         'SDC_IMAP_Lo/shared_code/pro_star/prostar'];
% CDF  folder
c_fld = 'matlab_cdf391_patch-arm64/';

addpath([s_fld d_fld]);
addpath([s_fld c_fld]);

%==============================================
%% Import & Process CDFs
close all; clear dir;

d_fld = ['/Users/mitchellshen/Library/CloudStorage/Dropbox-Personal/' ...
         'SDC_IMAP_Lo/shared_code/pro_star/prostar'];

f_fld = ['/Users/mitchellshen/Library/CloudStorage/Dropbox-Personal/' ...
         'SDC_IMAP_Lo/shared_code/pro_star/quicklook_MMS'];

k_fld = ['/Users/mitchellshen/Library/CloudStorage/Dropbox-Personal/' ...
         'SDC_IMAP_Lo/shared_code/pro_star/nhk'];

% Output folder for Figure 1 & 3 aves
outdir_f1 = fullfile(f_fld, 'fig_star_daily_f1');
if ~exist(outdir_f1, 'dir'); mkdir(outdir_f1); end

outdir_f3 = fullfile(f_fld, 'fig_star_daily_f3');
if ~exist(outdir_f3, 'dir'); mkdir(outdir_f3); end

files = dir(fullfile(d_fld, '*.cdf'));
names = {files.name};

% Parse: repoint##### and v###
tok = regexp(names, 'repoint(\d+)_v(\d+)\.cdf$', 'tokens', 'once');

% Keep only successfully-parsed ones
is_ok = ~cellfun(@isempty, tok);
files = files(is_ok);
tok   = tok(is_ok);

repoint = cellfun(@(t) str2double(t{1}), tok);
ver     = cellfun(@(t) str2double(t{2}), tok);

% For each repoint, choose max version (PROSTAR)
[uniq_rp, ~, grp] = unique(repoint);
best_idx = zeros(size(uniq_rp));

for k = 1:numel(uniq_rp)
    idx = find(grp == k);
    [~, imax] = max(ver(idx));
    best_idx(k) = idx(imax);
end

files_best = files(best_idx);

% Optional: sort by repoint (or by filename/date)
[~, srt] = sort(repoint(best_idx));
files_best = files_best(srt);

fprintf('Total CDFs: %d  ->  Highest-version-per-repoint: %d\n', ...
        numel(names), numel(files_best));

% Figure 2: initialize once (clean) + legend bookkeeping
figure(2); clf; hold on
label = {};                 % reset each run
h2 = gobjects(0);           % store line handles
j = 0;                      % compact index for labels/handles

% matrix for pivot angle
pivot = nan(1, numel(files_best));

for i = 1:numel(files_best)

    %% --- PROSTAR load ---
    fname = fullfile(d_fld, files_best(i).name);
    fprintf('Loading PROSTAR %s\n', files_best(i).name);

    [l1b_star_data, l1b_star_info] = spdfcdfread(fname, 'Structure', true);
    l1b_star_var = l1b_star_info.Variables(:, [1 4]); %#ok<NASGU>

    l1b_star_spin_ang_bin = (l1b_star_data(1).Data); %#ok<NASGU>
    l1b_star_met          = (l1b_star_data(2).Data); %#ok<NASGU>
    l1b_star_data_ave_amp = (l1b_star_data(3).Data);
    l1b_star_data_ave_amp(abs(l1b_star_data_ave_amp) > 1e30) = NaN; % fill values
    l1b_star_cunt_per_bin = (l1b_star_data(4).Data); %#ok<NASGU>
    l1b_star_epoch        = (l1b_star_data(5).Data);
    l1b_star_data_spn_ang = (l1b_star_data(6).Data);

    % Extract date / repoint / version from PROSTAR filename
    tok1 = regexp(files_best(i).name, '_(\d{8})-repoint(\d+)_v(\d+)', 'tokens', 'once');
    if ~isempty(tok1)
        date_str    = tok1{1};
        repoint_str = tok1{2};
        ver_str     = tok1{3};
    else
        tok_date = regexp(files_best(i).name, '_(\d{8})-', 'tokens', 'once');
        if ~isempty(tok_date), date_str = tok_date{1}; else, date_str = sprintf('idx%03d', i); end
        repoint_str = 'unknown';
        ver_str     = 'unknown';
    end

    %% --- NHK: find matching file (same date+repoint, any version) ---
    % Pattern: imap_lo_l1b_nhk_YYYYMMDD-repoint#####_v*.cdf
    nhk_pat  = sprintf('imap_lo_l1b_nhk_%s-repoint%s_v*.cdf', date_str, repoint_str);
    nhk_list = dir(fullfile(k_fld, nhk_pat));

    l1b_nhk_data = [];
    l1b_nhk_info = [];
    nhk_name = '';

    if isempty(nhk_list)
        warning('No NHK match for %s (pattern: %s)', files_best(i).name, nhk_pat);
    else
        % pick highest version among matches
        nhk_names = {nhk_list.name};
        tv = regexp(nhk_names, '_v(\d+)\.cdf$', 'tokens', 'once');
        v  = cellfun(@(x) str2double(x{1}), tv);
        [~, imax] = max(v);

        nhk_name = nhk_names{imax};
        kname = fullfile(k_fld, nhk_name);

        fprintf('Loading NHK    %s\n', nhk_name);
        [l1b_nhk_data, l1b_nhk_info] = spdfcdfread(kname, 'Structure', true);
    end
    
        % load pivot angle
        varname = 'pcc_cumulative_cnt_pri';
        idx = find(strcmp(l1b_nhk_info.Variables(:,1), varname), 1);  % find the row
        
        if isempty(idx)
            warning('NHK variable not found: %s in %s', varname, kname);
            pivot(i) = NaN;
        else
            l1b_nhk_data_pivot = l1b_nhk_data(idx).Data;
        
            % Optional: remove fill values if they are huge
            l1b_nhk_data_pivot(abs(l1b_nhk_data_pivot) > 1e30) = NaN;
            pivot(i) = round(mean(l1b_nhk_data_pivot, 'omitnan'));
        end


    %% Figure 1: Ave. STAR oplot (save per file)
    f1 = figure(1); clf; hold on;
    set(f1, 'Visible', 'off');

    plot(l1b_star_data_spn_ang, l1b_star_data_ave_amp);
    plot(l1b_star_data_spn_ang, mean(l1b_star_data_ave_amp), 'k-', 'LineWidth', 1.6);

    xlabel('NEP Angle (degree)');
    ylabel('Star Sensor Amplitude (mV)');
    xlim([0  360]); xticks([0:30:360]);
    ylim([0 2400]); yticks([0:200:2400]);

    labels1 = [arrayfun(@(k) sprintf('spin_ave_#%d',k), ...
        1:size(l1b_star_data_ave_amp,1), 'UniformOutput', false), ...
        {'Ave of L1b star'}];
    lgd1 = legend(labels1,'Location','eastoutside', 'Interpreter','none'); %#ok<NASGU>
    lgd1.NumColumns = 2;
    title(['Star Amplitude ' date_str ', Repoint ' repoint_str ', Pivot ' num2str(pivot(i),'%03d')])
    set(gca,'lineWidth',1,'FontSize',12,'FontWeight','default');
    grid on; box on;

    % Save Figure 1 with 300 dpi
    f1name = sprintf('F1_star_amp_%s_repoint%s_v%s.png', date_str, repoint_str, ver_str);
    exportgraphics(f1, fullfile(outdir_f1, f1name), 'Resolution', 300);

    %% Figure 2: plot mean curve per day (collect handles/labels for final legend)
    f2 = figure(2);
    j = j + 1;

    h2(j) = plot(l1b_star_data_spn_ang, mean(l1b_star_data_ave_amp), ...
                 'LineStyle','-', 'LineWidth', 1);

    % If you want to show NHK version in legend, you can add it here:
    % label{j} = sprintf('%s, pivot (nhk:%s)', date_str, nhk_name);
    label{j} = sprintf('%s, \x03B5 %03d', date_str, pivot(i));  % ε and 3-digit zero pad

    xlabel('NEP Angle (degree)');
    ylabel('Star Sensor Amplitude (mV)');
    xlim([0  360]); xticks([0:30:360]);
    ylim([0 2400]); yticks([0:200:2400]);

    set(gca,'lineWidth',1,'FontSize',12,'FontWeight','default');
    grid on; box on;

    %% Figure 3: 2D colormap (save per file)
    t = datetime(l1b_star_epoch,'ConvertFrom','datenum','TimeZone','UTC');
    Z = l1b_star_data_ave_amp.';   % (spin angle × time)

    f3 = figure(3);
    set(f3, 'Visible', 'off');
    h = imagesc(t, l1b_star_data_spn_ang, Z);
    axis xy
    ylim([0 360]); yticks([0:30:360]);
    hcb = colorbar('v');
    clim([0 2400]);
    set(hcb,'YTick',[0:200:2400])
    colormap("parula");
    xlabel('Epoch');
    ylabel('NEP Angle (degree)');
    ylabel(hcb,'Star Sensor Amplitude (mV)','FontSize',12);
    title(['Star Spectra ' date_str ', Repoint ' repoint_str ', Pivot ' num2str(pivot(i),'%03d')])


    set(h,'AlphaData', ~isnan(Z))
    set(gca,'Color','w');
    set(gca,'lineWidth',1,'FontSize',12,'FontWeight','default');
    box on;

    f3name = sprintf('F3_star_spc_%s_repoint%s_v%s.png', date_str, repoint_str, ver_str);
    exportgraphics(f3, fullfile(outdir_f3, f3name), 'Resolution', 300);

    %% ---- PLACEHOLDER: do your NHK plotting here ----
    % if ~isempty(l1b_nhk_data)
    %     % Example: inspect available variables
    %     % disp({l1b_nhk_data.Name})  % if it's struct-style; otherwise use l1b_nhk_info.Variables
    % end

end

% Final legend for Figure 2
f2 = figure(2);
lgd2 = legend(h2, label,'Location','eastoutside', 'Interpreter','none');
lgd2.NumColumns = 2;

f2name = sprintf('F2_star_combo.png');
    exportgraphics(f2, fullfile(f_fld, f2name), 'Resolution', 300);