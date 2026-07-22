%% MAIN
close all; clearvars; clc; format long; warning('off','all')
datetime.setDefaultFormats('default','yyyy-MM-dd hh:mm:ss')
addpath('~/Users/mitchellshen/mice/src/mice')
addpath('~/Users/mitchellshen/mice/lib')
savepath
rehash toolboxcache
tls_kernel = '/Users/mitchellshen/CAVA_IDL/kernels/lsk//naif0012.tls';
spk_kernel = '/Users/mitchellshen/mice//de440s.bsp';

%
pivot_list = [75 90 105];
species = 'H';

fld = '/Users/mitchellshen/Library/CloudStorage/Dropbox-Personal/';

dataset = struct();

for ipiv = 1:length(pivot_list)

    pivot = pivot_list(ipiv);

    fprintf('\nProcessing pivot %d ...\n', pivot);

    qmp = ['IMAP_ILO/SCI_ISN/p' num2str(pivot,'%03d') '_daily_' species];
    data_dir = fullfile(fld, qmp);

    % Get all CSV files
    files = dir(fullfile(data_dir, '*.csv'));

    if isempty(files)
        warning('No CSV files found for pivot %d in:\n%s', pivot, data_dir);
        continue
    end

    % Sort by filename
    [~, idx] = sort({files.name});
    files = files(idx);

    % Initialize output structure
    all_data = struct();

    % Import all csv for this pivot
    for ifile = 1:length(files)
        fname = files(ifile).name;
        fpath = fullfile(files(ifile).folder, fname);

        % Read CSV
        M = readmatrix(fpath);

        % Parse year-day and ESA step from filename
        tok = regexp(fname, 'data_YD_(\d+)_esa(\d+)\.csv', 'tokens');

        if isempty(tok)
            warning('Filename not matched: %s', fname);
            continue
        end

        yd  = tok{1}{1};
        esa = str2double(tok{1}{2});

        yd_field  = ['YD_' yd];
        esa_field = ['esa' num2str(esa)];

        all_data.(yd_field).(esa_field) = M;
    end

    % Flatten all to one array
    meta = [];

    yd_list = fieldnames(all_data);

    for iyd = 1:length(yd_list)
        yd_field = yd_list{iyd};
        yd_str = yd_field(4:end);
        yd_num = str2double(yd_str);

        esa_list = fieldnames(all_data.(yd_field));

        for iesa = 1:length(esa_list)
            esa_field = esa_list{iesa};
            esa_num = str2double(esa_field(4:end));

            M = all_data.(yd_field).(esa_field);
            N = size(M,1);

            yd_col  = yd_num * ones(N,1);
            esa_col = esa_num * ones(N,1);

            meta = [meta; M, yd_col, esa_col];
        end
    end

    % Average first 30 rows for each (yd, esa)
    col_bins     = 1;
    col_counts   = 2;
    col_ra       = 3;
    col_dec      = 4;
    col_expo     = 5;
    col_spin_ra  = 6;
    col_spin_dec = 7;
    col_yd       = 8;
    col_esa      = 9;

    yd_vals  = unique(meta(:,col_yd));
    esa_vals = unique(meta(:,col_esa));

    ram = [];

    for iyd = 1:length(yd_vals)
        yd_now = yd_vals(iyd);

        for iesa = 1:length(esa_vals)
            esa_now = esa_vals(iesa);

            idx = meta(:,col_yd) == yd_now & meta(:,col_esa) == esa_now;
            Msub = meta(idx,:);

            if isempty(Msub)
                continue
            end

            n_take = min(30, size(Msub,1));
            M30 = Msub(1:n_take,:);
            % M30 = Msub(11:21,:); % close to ecliptic
            % M30 = Msub(1:60,:); % both ram + aram


            counts_avg   = mean(M30(:,col_counts),   'omitnan');
            dec_avg      = mean(M30(:,col_dec),      'omitnan');
            expo_avg     = mean(M30(:,col_expo),     'omitnan');
            spin_dec_avg = mean(M30(:,col_spin_dec), 'omitnan');
            rate_ave     = sum(M30(:,col_counts), 'omitnan') ./ ...
                           sum(M30(:,col_expo),   'omitnan');

            ra_avg      = circ_mean_deg(M30(:,col_ra));
            spin_ra_avg = circ_mean_deg(M30(:,col_spin_ra));

            bin_start = 1;
            bin_end   = n_take;

            ram = [ram;
                   yd_now, esa_now, bin_start, bin_end, ...
                   counts_avg, ra_avg, dec_avg, expo_avg, ...
                   spin_ra_avg, spin_dec_avg, rate_ave];
        end
    end

    % Convert yyyydoy to datetime
    yd = ram(:,1);
    year = floor(yd / 1000);
    doy  = mod(yd, 1000);

    time = datetime(year,1,1) + days(doy - 1);

    % Store everything for this pivot
    dataset(ipiv).pivot    = pivot;
    dataset(ipiv).species  = species;
    dataset(ipiv).data_dir = data_dir;
    dataset(ipiv).all_data = all_data;
    dataset(ipiv).meta     = meta;
    dataset(ipiv).ram      = ram;
    dataset(ipiv).time     = time;

    fprintf('Done pivot %d: %d averaged rows\n', pivot, size(ram,1));
end

%% PLOT

figure(1); clf
hold on

esa_plot = [1:7];

for k = 1:length(esa_plot)
    esa_now = esa_plot(k);

    idx = ram(:,2) == esa_now;

    plot(time(idx), ram(idx,11), 'o', 'LineWidth', 1.5, ...
        'DisplayName', ['ESA ' num2str(esa_now)]);
end

xlabel('Time')
ylabel('Rate average')
title('Rate average vs time')
legend('Location','best')
grid on
box on




%% ORBIT PLOT
% dataset(ipiv).pivot
% dataset(ipiv).ram
% dataset(ipiv).time
%
% ram columns:
% 1  yyyydoy
% 2  esa
% 3  bin_start
% 4  bin_end
% 5  counts_avg
% 6  ra_avg
% 7  dec_avg
% 8  expo_avg
% 9  spin_ra_avg
% 10 spin_dec_avg
% 11 rate_ave

figure(3); clf
hold on

% Upstream
% Put interstellar He upstream on +X axis
lon_upstream_deg = 75.4 +180;
lon_upstream = deg2rad(lon_upstream_deg);


% Reference epoch for orbital angle
% Sep 22 --> 3 o'clock
% t_ref = datetime(2025,6,1,0,0,0,0);
% t_ref = datetime(2025,6,9,0,0,0,0);
% t_ref = datetime(2025,9,23,0,0,0,0);


% Pivot groups: one center radius per pivot
% Example:
%   pivot 75  -> inner group
%   pivot 90  -> middle group
%   pivot 105 -> outer group
pivot_centers = [0.75 1.00 1.25];

% Small ESA offsets around each pivot center
esa_offset = [-0.09 -0.06 -0.03 0 0.03 0.06 0.09];

% Orbit circle angle
th = linspace(0,2*pi,400);

% Draw reference circles for all pivots and ESA steps
for ipiv = 1:length(dataset)

    if ~isfield(dataset(ipiv),'ram') || isempty(dataset(ipiv).ram)
        continue
    end

    r_center = pivot_centers(ipiv);
    esa_radii = r_center + esa_offset;

    for esa_now = 1:7
        r0 = esa_radii(esa_now);

        % Highlight ESA 4 ring for each pivot
        if esa_now == 4
            plot(r0*cos(th), r0*sin(th), 'k-', 'LineWidth', 1.0, ...
                'HandleVisibility','off');
        else
            plot(r0*cos(th), r0*sin(th), '-', ...
                'Color', [0.82 0.82 0.82], ...
                'LineWidth', 0.7, ...
                'HandleVisibility','off');
        end
    end
end

% Month markers on the pivot=90, ESA=4 reference circle
% Here I assume pivot_centers(2) corresponds to pivot 90
% month_dates = t_ref + calmonths(0:11);
% theta_m = 2*pi * days(month_dates - t_ref) / 365.25;

month_dates = [ ...
    datetime(2025,7,1) + calmonths(0:6), ...
    datetime(2026,1,1) + calmonths(0:5) ]';
lon_m = get_earth_lon_eclipj2000(month_dates, tls_kernel, spk_kernel);
lon_m_plot = mod(lon_m - lon_upstream, 2*pi);

% Month labels outside the outermost pivot/ESA ring
label_r = max(pivot_centers + max(esa_offset)) + 0.18;
for k = 1:length(lon_m)
    text(label_r*cos(lon_m_plot(k)), label_r*sin(lon_m_plot(k)), ...
        datestr(month_dates(k),'mmm'), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'FontSize',10);
end

% Instead of doing dots, we then do thin lines
% r_month = pivot_centers(2) + 0;   % pivot 90, ESA 4
% plot(r_month*cos(theta_m), r_month*sin(theta_m), 'm.', ...
%     'MarkerSize', 12, 'HandleVisibility','off')

for k = 1:length(lon_m)

    % default style
    lw = 0.6;
    col = [0.7 0.7 0.7];

    label = datestr(month_dates(k),'mmm');

    % if strcmp(label,'Sep')
    %     lw = 1.5;
    %     col = [0 0 0];   % bold reference axis
    % elseif strcmp(label,'Dec')
    %     lw = 1.5;
    %     col = [1 0 0];   % highlight cone direction
    % end

    % plot([0, (label_r-0.05)*cos(theta_m(k))], ...
    %      [0, (label_r-0.05)*sin(theta_m(k))], ...
    plot([0, (label_r-0.05)*cos(lon_m_plot(k))], ...
         [0, (label_r-0.05)*sin(lon_m_plot(k))], ...)
         '-', 'Color', col, 'LineWidth', lw, ...
         'HandleVisibility','off');
end

% Plot data for each pivot and ESA
for ipiv = 1:length(dataset)

    if ~isfield(dataset(ipiv),'ram') || isempty(dataset(ipiv).ram)
        continue
    end

    pivot_now = dataset(ipiv).pivot;
    ram  = dataset(ipiv).ram;
    time = dataset(ipiv).time;

    r_center = pivot_centers(ipiv);
    esa_radii = r_center + esa_offset;

    for esa_now = 1:7
        idx = ram(:,2) == esa_now;

        if ~any(idx)
            continue
        end

        time_sel = time(idx);
        rate_sel = ram(idx,11);

        % Sort by time
        [time_sel, isrt] = sort(time_sel);
        rate_sel = rate_sel(isrt);

        % Orbital phase
        % theta = 2*pi * days(time_sel - t_ref) / 365.25;
        lon_sel = get_earth_lon_eclipj2000(time_sel, tls_kernel, spk_kernel);
        lon_sel_plot = mod(lon_sel - lon_upstream, 2*pi);

        % Radius for this pivot + ESA
        r0 = esa_radii(esa_now);

        % x = r0 * cos(theta);
        % y = r0 * sin(theta);
        x = r0 * cos(lon_sel_plot);
        y = r0 * sin(lon_sel_plot);


        scatter(x, y, 15, rate_sel, 'filled', ...
            'DisplayName', sprintf('P%d ESA%d', pivot_now, esa_now));
    end
end

% Sun
plot(0,0,'yo','MarkerFaceColor',[1 0.5 0], 'MarkerEdgeColor',[1 0.5 0], ...
    'MarkerSize',16, ...
    'DisplayName','Sun')

% Pivot annotation
for ipiv = 1:length(pivot_list)
    str(ipiv,:) = ['p' num2str(pivot_list(ipiv),'%03d')];
end
text(zeros(1,size(str,1)),pivot_centers-0.04,str,...
    'HorizontalAlignment','center','Color','b');
% text(pivot_centers-0.04,zeros(1,size(str,1)),str,...
%     'HorizontalAlignment','center','Color','b');


axis equal
xlim(1.15*[-label_r label_r])
ylim(1.15*[-label_r label_r])
grid on
box on
xlabel('X (AU)')
ylabel('Y (AU)')
title('Circular orbit view by pivot and ESA')
cb = colorbar;
cb.Label.String = 'Rate average';
colormap("turbo")
% set(gca, 'ColorScale', 'log')

% legend('Location','eastoutside')



% Overlay SWAPI PUI on inner ring

pui_path = [fld 'IMAP_ILO/SCI_ISN/swapi_pui_density.mat'];
S = load(pui_path);

% If fields are top-level in the MAT file:
if isfield(S,'t_mid_arr')
    time_pui      = S.t_mid_arr;
    dens_pui      = double(S.pui_dens_arr);
    dens_pui_norm = double(S.pui_dens_arr_norm);

% If MAT file contains a struct named PUI:
elseif isfield(S,'PUI')
    time_pui      = S.PUI.t_mid_arr;
    dens_pui      = double(S.PUI.pui_dens_arr);
    dens_pui_norm = double(S.PUI.pui_dens_arr_norm);

else
    error('Could not find PUI variables in swapi_pui_density.mat')
end

% Keep only finite values
good = isfinite(dens_pui_norm) & ~isnat(time_pui);
time_pui      = time_pui(good);
dens_pui      = dens_pui(good);
dens_pui_norm = dens_pui_norm(good);

% Convert to same orbital angle frame as IMAP-Lo
% t_ref = datetime(t_ref,'TimeZone',time_pui.TimeZone);
% theta_pui = 2*pi * days(time_pui - t_ref) / 365.25;
% theta_pui = mod(theta_pui, 2*pi);
lon_pui = get_earth_lon_eclipj2000(time_pui, tls_kernel, spk_kernel);
lon_pui_plot = mod(lon_pui - lon_upstream, 2*pi);

% Sort by angle for a smooth trace
[lon_pui, isrt] = sort(lon_pui);
dens_pui_norm = dens_pui_norm(isrt);

% --- Choose the PUI reference ring radius ---
r_pui0 = 0.3;   % adjust if you want it slightly larger/smaller

% Draw faint PUI reference circle
th = linspace(0,2*pi,400);
plot(r_pui0*cos(th), r_pui0*sin(th), '-', ...
    'Color', [0.85 0.65 0.90], ...
    'LineWidth', 0.5, ...
    'HandleVisibility','off');

% --- Convert normalized density into a small radial modulation ---
% rescale to [0,1]
dmin = min(dens_pui_norm);
dmax = max(dens_pui_norm);

if dmax > dmin
    dens01 = (dens_pui_norm - dmin) ./ (dmax - dmin);
else
    dens01 = zeros(size(dens_pui_norm));
end

% radial excursion size
dr_amp = 0.8;   % try 0.06 to 0.12 depending on how strong you want it

r_pui = r_pui0 + dr_amp * dens01;

% x_pui = r_pui .* cos(theta_pui);
% y_pui = r_pui .* sin(theta_pui);
x_pui = r_pui .* cos(lon_pui_plot);
y_pui = r_pui .* sin(lon_pui_plot);

% Plot the thin magenta PUI curve
% plot(x_pui, y_pui, '-', ...
%     'Color', [0.95 0.25 0.95], ...
%     'LineWidth', 1.6, ...
%     'DisplayName', 'SWAPI PUI norm');
scatter(x_pui, y_pui, 1, 'filled', ...
    'MarkerFaceColor', [0.95 0.25 0.95], ...
    'DisplayName', 'SWAPI PUI norm');

set(gca,'FontSize',12,'FontWeight','default');



%% ORBIT PLOT2
% dataset(ipiv).pivot
% dataset(ipiv).ram
% dataset(ipiv).time
%
% ram columns:
% 1  yyyydoy
% 2  esa
% 3  bin_start
% 4  bin_end
% 5  counts_avg
% 6  ra_avg
% 7  dec_avg
% 8  expo_avg
% 9  spin_ra_avg
% 10 spin_dec_avg
% 11 rate_ave

figure(4); clf
hold on

% Upstream
% Put interstellar He upstream on +X axis
lon_upstream_deg = 75.4 + 180
lon_upstream = deg2rad(lon_upstream_deg);

% Two IMAP-Lo rings only:
% inner ring = combined p075 + p105
% outer ring = p090
pivot_centers = [0.85 1.10];   % [inner, outer]
esa_now = 4;                   % ESA 4 only

% Orbit circle angle
th = linspace(0,2*pi,400);

% Draw only 2 reference circles for IMAP-Lo
for ir = 1:2
    r0 = pivot_centers(ir);
    plot(r0*cos(th), r0*sin(th), 'k-', 'LineWidth', 1.0, ...
        'HandleVisibility','off');
end

% Month annotations
month_dates = [ ...
    datetime(2025,7,1) + calmonths(0:6), ...
    datetime(2026,1,1) + calmonths(0:5) ]';
lon_m = get_earth_lon_eclipj2000(month_dates, tls_kernel, spk_kernel);
lon_m_plot = mod(lon_m - lon_upstream, 2*pi);

% Month labels outside the outermost ring
label_r = max(pivot_centers) + 0.18;
for k = 1:length(lon_m)
    text(label_r*cos(lon_m_plot(k)), label_r*sin(lon_m_plot(k)), ...
        datestr(month_dates(k),'mmm'), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'FontSize',10);
end

% Thin month spokes
for k = 1:length(lon_m)

    lw = 0.6;
    col = [0.7 0.7 0.7];

    plot([0, (label_r-0.05)*cos(lon_m_plot(k))], ...
         [0, (label_r-0.05)*sin(lon_m_plot(k))], ...
         '-', 'Color', col, 'LineWidth', lw, ...
         'HandleVisibility','off');
end

% Plot IMAP-Lo data: ESA 4 only
for ipiv = 1:length(dataset)

    if ~isfield(dataset(ipiv),'ram') || isempty(dataset(ipiv).ram)
        continue
    end

    pivot_now = dataset(ipiv).pivot;
    ram  = dataset(ipiv).ram;
    time = dataset(ipiv).time;

    idx = ram(:,2) == esa_now;

    if ~any(idx)
        continue
    end

    time_sel = time(idx);
    rate_sel = ram(idx,11);

    % Sort by time
    [time_sel, isrt] = sort(time_sel);
    rate_sel = rate_sel(isrt);

    lon_sel = get_earth_lon_eclipj2000(time_sel, tls_kernel, spk_kernel);
    lon_sel_plot = mod(lon_sel - lon_upstream, 2*pi);

    % Ring assignment
    if pivot_now == 90
        r0 = pivot_centers(2);         % outer ring
    elseif pivot_now == 75 || pivot_now == 105
        r0 = pivot_centers(1);         % inner ring
    else
        continue
    end

    x = r0 * cos(lon_sel_plot);
    y = r0 * sin(lon_sel_plot);

    scatter(x, y, 30, rate_sel, 'filled', ...
        'DisplayName', sprintf('P%d ESA4', pivot_now));
end

% Sun
plot(0,0,'yo','MarkerFaceColor',[1 0.5 0], 'MarkerEdgeColor',[1 0.5 0], ...
    'MarkerSize',16, ...
    'DisplayName','Sun')

% Pivot annotation
text(0, pivot_centers(1)-0.04, 'p075 + p105', ...
    'HorizontalAlignment','center','Color','b');

text(0, pivot_centers(2)-0.04, 'p090', ...
    'HorizontalAlignment','center','Color','b');

axis equal
xlim(1.15*[-label_r label_r])
ylim(1.15*[-label_r label_r])
grid on
box on
xlabel('X (AU)')
ylabel('Y (AU)')
title('Circular orbit view by pivot (ESA 4 only)')
cb = colorbar;
cb.Label.String = 'Rate average';
colormap("turbo")
% set(gca, 'ColorScale', 'log')

% legend('Location','eastoutside')

% Overlay SWAPI PUI on inner ring
% Keep SWAPI PUI logic the same

pui_path = [fld 'IMAP_ILO/SCI_ISN/swapi_pui_density.mat'];
S = load(pui_path);

% If fields are top-level in the MAT file:
if isfield(S,'t_mid_arr')
    time_pui      = S.t_mid_arr;
    dens_pui      = double(S.pui_dens_arr);
    dens_pui_norm = double(S.pui_dens_arr_norm);

% If MAT file contains a struct named PUI:
elseif isfield(S,'PUI')
    time_pui      = S.PUI.t_mid_arr;
    dens_pui      = double(S.PUI.pui_dens_arr);
    dens_pui_norm = double(S.PUI.pui_dens_arr_norm);

else
    error('Could not find PUI variables in swapi_pui_density.mat')
end

% Keep only finite values
good = isfinite(dens_pui_norm) & ~isnat(time_pui);
time_pui      = time_pui(good);
dens_pui      = dens_pui(good);
dens_pui_norm = dens_pui_norm(good);

% Convert to same orbital angle frame as IMAP-Lo
lon_pui = get_earth_lon_eclipj2000(time_pui, tls_kernel, spk_kernel);
lon_pui_plot = mod(lon_pui - lon_upstream, 2*pi);

% Sort by angle for a smooth trace
[lon_pui_plot, isrt] = sort(lon_pui_plot);
dens_pui_norm = dens_pui_norm(isrt);
dens_pui      = dens_pui(isrt);
time_pui      = time_pui(isrt);

% --- Choose the PUI reference ring radius ---
r_pui0 = 0.3;   % adjust if you want it slightly larger/smaller

% Draw faint PUI reference circle
th = linspace(0,2*pi,400);
plot(r_pui0*cos(th), r_pui0*sin(th), '-', ...
    'Color', [0.85 0.65 0.90], ...
    'LineWidth', 0.5, ...
    'HandleVisibility','off');

% --- Convert normalized density into a small radial modulation ---
% rescale to [0,1]
dmin = min(dens_pui_norm);
dmax = max(dens_pui_norm);

if dmax > dmin
    dens01 = (dens_pui_norm - dmin) ./ (dmax - dmin);
else
    dens01 = zeros(size(dens_pui_norm));
end

% radial excursion size
dr_amp = 0.8;   % keep your current setting

r_pui = r_pui0 + dr_amp * dens01;

x_pui = r_pui .* cos(lon_pui_plot);
y_pui = r_pui .* sin(lon_pui_plot);

% Plot the thin magenta PUI curve
% plot(x_pui, y_pui, '-', ...
%     'Color', [0.95 0.25 0.95], ...
%     'LineWidth', 1.6, ...
%     'DisplayName', 'SWAPI PUI norm');
scatter(x_pui, y_pui, 1, 'filled', ...
    'MarkerFaceColor', [0.95 0.25 0.95], ...
    'DisplayName', 'SWAPI PUI norm');

set(gca,'FontSize',12,'FontWeight','default');

%% FIT AND ANNOTATE PEAKS
% Uses ESA 4 only, same upstream-rotated longitude convention as the plot

figure(4); hold on

% ---------- Gather ESA4 data ----------
theta_p90 = [];
rate_p90  = [];

theta_in  = [];
rate_in   = [];

for ipiv = 1:length(dataset)

    if ~isfield(dataset(ipiv),'ram') || isempty(dataset(ipiv).ram)
        continue
    end

    pivot_now = dataset(ipiv).pivot;
    ram  = dataset(ipiv).ram;
    time = dataset(ipiv).time;

    idx = ram(:,2) == 4;   % ESA 4 only
    if ~any(idx)
        continue
    end

    time_sel = time(idx);
    rate_sel = ram(idx,11);

    [time_sel, isrt] = sort(time_sel);
    rate_sel = rate_sel(isrt);

    lon_sel = get_earth_lon_eclipj2000(time_sel, tls_kernel, spk_kernel);
    lon_sel_plot = mod(lon_sel - lon_upstream, 2*pi);

    good = isfinite(lon_sel_plot) & isfinite(rate_sel);
    lon_sel_plot = lon_sel_plot(good);
    rate_sel = rate_sel(good);

    if pivot_now == 90
        theta_p90 = [theta_p90; lon_sel_plot(:)];
        rate_p90  = [rate_p90;  rate_sel(:)];
    elseif pivot_now == 75 || pivot_now == 105
        theta_in = [theta_in; lon_sel_plot(:)];
        rate_in  = [rate_in;  rate_sel(:)];
    end
end

% ---------- Sort for plotting ----------
[theta_p90, isrt] = sort(theta_p90);
rate_p90 = rate_p90(isrt);

[theta_in, isrt] = sort(theta_in);
rate_in = rate_in(isrt);

% ---------- Fit p090 with one wrapped Gaussian ----------
b0  = min(rate_p90);
[a0, imx] = max(rate_p90);
a0  = a0 - b0;
mu0 = theta_p90(imx);
s0  = deg2rad(18);

p0_1g = [b0, a0, mu0, s0];

opts = optimset('Display','off','MaxFunEvals',8000,'MaxIter',8000);
pfit_1g = fminsearch(@(p) chisq_1g(p, theta_p90, rate_p90), p0_1g, opts);

b_90   = pfit_1g(1);
a_90   = pfit_1g(2);
mu_90  = mod(pfit_1g(3), 2*pi);
sig_90 = abs(pfit_1g(4));

% ---------- Fit inner ring with two wrapped Gaussians ----------
% Crude initial guesses from the two strongest points separated by angle
[~, isrt_pk] = sort(rate_in, 'descend');
mu1_0 = theta_in(isrt_pk(1));

mu2_0 = [];
for ii = 2:length(isrt_pk)
    cand = theta_in(isrt_pk(ii));
    if abs(wrapdiff(cand, mu1_0)) > deg2rad(15)
        mu2_0 = cand;
        break
    end
end
if isempty(mu2_0)
    mu2_0 = mod(mu1_0 + deg2rad(25), 2*pi);
end

b0_in = min(rate_in);
a1_0  = max(rate_in) - b0_in;
a2_0  = 0.7 * a1_0;
s1_0  = deg2rad(12);
s2_0  = deg2rad(12);

p0_2g = [b0_in, a1_0, mu1_0, s1_0, a2_0, mu2_0, s2_0];

pfit_2g = fminsearch(@(p) chisq_2g(p, theta_in, rate_in), p0_2g, opts);

b_in   = pfit_2g(1);
a1_in  = pfit_2g(2);
mu1_in = mod(pfit_2g(3), 2*pi);
s1_in  = abs(pfit_2g(4));
a2_in  = pfit_2g(5);
mu2_in = mod(pfit_2g(6), 2*pi);
s2_in  = abs(pfit_2g(7));

% Sort the two inner peaks by angle for cleaner labeling
mu_pair = [mu1_in, mu2_in];
a_pair  = [a1_in,  a2_in];
s_pair  = [s1_in,  s2_in];

[mu_pair, ord] = sort(mu_pair);
a_pair = a_pair(ord);
s_pair = s_pair(ord);

mu1_in = mu_pair(1);  mu2_in = mu_pair(2);
a1_in  = a_pair(1);   a2_in  = a_pair(2);
s1_in  = s_pair(1);   s2_in  = s_pair(2);

% ---------- Annotate peaks on orbit plot ----------
r_inner = pivot_centers(1);
r_outer = pivot_centers(2);

% p090 single peak
plot([0, r_outer*cos(mu_90)], [0, r_outer*sin(mu_90)], ...
    '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.2, 'HandleVisibility','off');
plot(r_outer*cos(mu_90), r_outer*sin(mu_90), 'o', ...
    'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'k', ...
    'MarkerSize', 8, 'HandleVisibility','off');

text((r_outer+0.08)*cos(mu_90), (r_outer+0.08)*sin(mu_90), ...
    sprintf('%.1fº', rad2deg(mu_90)), ...
    'HorizontalAlignment','center', 'Color','k', 'FontSize',10);

% inner ring two peaks
mu_in_all = [mu1_in, mu2_in];
for k = 1:2
    mu_now = mu_in_all(k);

    plot([0, r_inner*cos(mu_now)], [0, r_inner*sin(mu_now)], ...
        '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.0, ...
        'HandleVisibility','off');

    plot(r_inner*cos(mu_now), r_inner*sin(mu_now), 's', ...
        'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'k', ...
        'MarkerSize', 7, 'HandleVisibility','off');

    text((r_inner+0.08)*cos(mu_now), (r_inner+0.08)*sin(mu_now), ...
        sprintf('%.1fº', rad2deg(mu_now)), ...
        'HorizontalAlignment','center', 'Color','k', 'FontSize',9);
end

% ---------- Optional diagnostic figure ----------
theta_grid = linspace(0, 2*pi, 1000).';

y90_fit = model_1g(pfit_1g, theta_grid);
yin_fit = model_2g(pfit_2g, theta_grid);

figure(5); clf

subplot(2,1,1)
plot(rad2deg(theta_p90), rate_p90, 'ko', 'MarkerFaceColor', [0.3 0.6 1]); hold on
plot(rad2deg(theta_grid), y90_fit, 'k-', 'LineWidth', 2)
xline(rad2deg(mu_90), '--c', 'LineWidth', 1.2)
xlim([0 360]); grid on; box on
ylabel('Rate average')
title('p090, ESA 4: one wrapped Gaussian fit')

subplot(2,1,2)
plot(rad2deg(theta_in), rate_in, 'ko', 'MarkerFaceColor', [1 0.4 1]); hold on
plot(rad2deg(theta_grid), yin_fit, 'k-', 'LineWidth', 2)
xline(rad2deg(mu1_in), '--m', 'LineWidth', 1.2)
xline(rad2deg(mu2_in), '--m', 'LineWidth', 1.2)
xlim([0 360]); grid on; box on
xlabel('Rotated orbital longitude (deg)')
ylabel('Rate average')
title('p075 + p105, ESA 4: two wrapped Gaussian fit')

%% ORBIT PLOT: ONE CIRCLE, SHOWING PEAKS DETERMINED FROM THE 2-RING FIT

figure(6); clf
hold on

% Upstream
lon_upstream_deg = 75.4 +180;
lon_upstream = deg2rad(lon_upstream_deg);

% Single IMAP-Lo ring
r_imap = 1.00;
th = linspace(0,2*pi,400);

plot(r_imap*cos(th), r_imap*sin(th), 'k-', 'LineWidth', 1.0, ...
    'HandleVisibility','off');

% Month annotations (same convention as before)
month_dates = [ ...
    datetime(2025,7,1) + calmonths(0:6), ...
    datetime(2026,1,1) + calmonths(0:5) ]';
lon_m = get_earth_lon_eclipj2000(month_dates, tls_kernel, spk_kernel);
lon_m_plot = mod(lon_m - lon_upstream, 2*pi);

label_r = r_imap + 0.22;

for k = 1:length(lon_m)
    text(label_r*cos(lon_m_plot(k)), label_r*sin(lon_m_plot(k)), ...
        datestr(month_dates(k),'mmm'), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'FontSize',10);
end

for k = 1:length(lon_m)
    plot([0, (label_r-0.05)*cos(lon_m_plot(k))], ...
         [0, (label_r-0.05)*sin(lon_m_plot(k))], ...
         '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.6, ...
         'HandleVisibility','off');
end

% Plot all 3 pivots together on one circle, ESA 4 only
esa_now = 4;

for ipiv = 1:length(dataset)

    if ~isfield(dataset(ipiv),'ram') || isempty(dataset(ipiv).ram)
        continue
    end

    ram  = dataset(ipiv).ram;
    time = dataset(ipiv).time;

    idx = ram(:,2) == esa_now;
    if ~any(idx)
        continue
    end

    time_sel = time(idx);
    rate_sel = ram(idx,11);

    [time_sel, isrt] = sort(time_sel);
    rate_sel = rate_sel(isrt);

    lon_sel = get_earth_lon_eclipj2000(time_sel, tls_kernel, spk_kernel);
    lon_sel_plot = mod(lon_sel - lon_upstream, 2*pi);

    x = r_imap * cos(lon_sel_plot);
    y = r_imap * sin(lon_sel_plot);

    scatter(x, y, 18, rate_sel, 'filled', ...
        'HandleVisibility','off');
end

% Sun
plot(0,0,'yo','MarkerFaceColor',[1 0.5 0], 'MarkerEdgeColor',[1 0.5 0], ...
    'MarkerSize',16, ...
    'DisplayName','Sun')

% Ring annotation
% text(0, r_imap-0.04, 'p075 + p090 + p105', ...
%     'HorizontalAlignment','center', 'Color','b');

axis equal
xlim(1.15*[-label_r label_r])
ylim(1.15*[-label_r label_r])
grid on
box on
xlabel('X (AU)')
ylabel('Y (AU)')
title('Circular orbit view by pivot (ESA 4 only, one circle)')
cb = colorbar;
cb.Label.String = 'Rate average';
colormap("turbo")

% Overlay SWAPI PUI on inner ring
% Keep your current SWAPI PUI handling unchanged

pui_path = [fld 'IMAP_ILO/SCI_ISN/swapi_pui_density.mat'];
S = load(pui_path);

if isfield(S,'t_mid_arr')
    time_pui      = S.t_mid_arr;
    dens_pui      = double(S.pui_dens_arr);
    dens_pui_norm = double(S.pui_dens_arr_norm);
elseif isfield(S,'PUI')
    time_pui      = S.PUI.t_mid_arr;
    dens_pui      = double(S.PUI.pui_dens_arr);
    dens_pui_norm = double(S.PUI.pui_dens_arr_norm);
else
    error('Could not find PUI variables in swapi_pui_density.mat')
end

good = isfinite(dens_pui_norm) & ~isnat(time_pui);
time_pui      = time_pui(good);
dens_pui      = dens_pui(good);
dens_pui_norm = dens_pui_norm(good);

lon_pui = get_earth_lon_eclipj2000(time_pui, tls_kernel, spk_kernel);
lon_pui_plot = mod(lon_pui - lon_upstream, 2*pi);

[lon_pui_plot, isrt] = sort(lon_pui_plot);
dens_pui_norm = dens_pui_norm(isrt);

r_pui0 = 0.3;

plot(r_pui0*cos(th), r_pui0*sin(th), '-', ...
    'Color', [0.85 0.65 0.90], ...
    'LineWidth', 0.5, ...
    'HandleVisibility','off');

dmin = min(dens_pui_norm);
dmax = max(dens_pui_norm);

if dmax > dmin
    dens01 = (dens_pui_norm - dmin) ./ (dmax - dmin);
else
    dens01 = zeros(size(dens_pui_norm));
end

dr_amp = 0.8;   % unchanged from your current SWAPI setting
r_pui = r_pui0 + dr_amp * dens01;

x_pui = r_pui .* cos(lon_pui_plot);
y_pui = r_pui .* sin(lon_pui_plot);

scatter(x_pui, y_pui, 1, 'filled', ...
    'MarkerFaceColor', [0.95 0.25 0.95], ...
    'DisplayName', 'SWAPI PUI norm');

% Annotate the peaks determined from the previous 2-ring plot
% Requires mu_90, mu1_in, mu2_in already existing

mu_show = [mu1_in, mu_90, mu2_in];
lab_show = {'p075','p090','p105'};
mk_show  = {'s','o','s'};
fc_show  = {'m','c','m'};
ls_show  = {'--','-','--'};

for k = 1:3
    mu_now = mu_show(k);

    % radial guide
    plot([0, r_imap*cos(mu_now)], [0, r_imap*sin(mu_now)], ...
        ls_show{k}, 'Color', [0.35 0.35 0.35], 'LineWidth', 1.0, ...
        'HandleVisibility','off');
  
    % marker at peak
    plot(r_imap*cos(mu_now), r_imap*sin(mu_now), mk_show{k}, ...
        'MarkerFaceColor', fc_show{k}, 'MarkerEdgeColor', 'k', ...
        'MarkerSize', 8, 'HandleVisibility','off');

    % label placement
    x_txt = (r_imap + 0.18) * cos(mu_now);
    y_txt = (r_imap + 0.18) * sin(mu_now);

    if cos(mu_now) > 0.2
        halign = 'left';
    elseif cos(mu_now) < -0.2
        halign = 'right';
    else
        halign = 'center';
    end

    text(x_txt, y_txt, ...
        sprintf('%s\n%.1f°', lab_show{k}, rad2deg(mu_now)), ...
        'HorizontalAlignment', halign, ...
        'VerticalAlignment', 'middle', ...
        'FontSize', 9, ...
        'BackgroundColor', 'w', ...
        'Margin', 2);
end

set(gca,'FontSize',12,'FontWeight','default');



%% SWAPI PUI test
figure(10)
subplot(2,1,1)
scatter(time_pui,dens_pui,5);
grid on; box on;
subplot(2,1,2)
scatter(time_pui,dens_pui_norm,5);
grid on; box on;



%% ===== BREAK =====


%% ESTIMATE CONE CENTER USING ALL PIVOTS TOGETHER

esa_fit = 3;

theta_all_fit = [];
rate_all_fit  = [];
pivot_id_fit  = [];

for ipiv = 1:length(dataset)

    if ~isfield(dataset(ipiv),'ram') || isempty(dataset(ipiv).ram)
        continue
    end

    ram_i  = dataset(ipiv).ram;
    time_i = dataset(ipiv).time;

    idx_fit = ram_i(:,2) == esa_fit;
    if ~any(idx_fit)
        continue
    end

    time_fit_i = time_i(idx_fit);
    rate_fit_i = ram_i(idx_fit,11);

    theta_fit_i = 2*pi * days(time_fit_i - t_ref) / 365.25;
    theta_fit_i = mod(theta_fit_i, 2*pi);

    good = isfinite(theta_fit_i) & isfinite(rate_fit_i);
    theta_fit_i = theta_fit_i(good);
    rate_fit_i  = rate_fit_i(good);

    theta_all_fit = [theta_all_fit; theta_fit_i];
    rate_all_fit  = [rate_all_fit;  rate_fit_i];
    pivot_id_fit  = [pivot_id_fit;  dataset(ipiv).pivot * ones(sum(good),1)];
end

% sort for cleaner plotting later
[theta_all_fit, isrt] = sort(theta_all_fit);
rate_all_fit = rate_all_fit(isrt);
pivot_id_fit = pivot_id_fit(isrt);

% Quick weighted center estimate
pct_use = 70;   % top 30% brightest points
thr = prctile(rate_all_fit, pct_use);
use = rate_all_fit >= thr;

theta_use = theta_all_fit(use);
w_use = rate_all_fit(use);

xw = sum(w_use .* cos(theta_use));
yw = sum(w_use .* sin(theta_use));

theta_w = atan2(yw, xw);
if theta_w < 0
    theta_w = theta_w + 2*pi;
end

% Gaussian-like fit in wrapped angle
bkg0   = min(rate_all_fit);
pmax   = max(rate_all_fit);
A0     = pmax - bkg0;
theta0 = theta_w;
sigma0 = deg2rad(20);

p0 = [bkg0, A0, theta0, sigma0];

opts = optimset('Display','off', 'MaxFunEvals',5000, 'MaxIter',5000);
pbest = fminsearch(@(p) cone_chisq(p, theta_all_fit, rate_all_fit), p0, opts);

bkg_best   = pbest(1);
A_best     = pbest(2);
theta_best = mod(pbest(3), 2*pi);
sigma_best = abs(pbest(4));
fwhm_best  = 2.355 * sigma_best;

fprintf('\n=== Combined cone estimate for ESA %d (all pivots) ===\n', esa_fit);
fprintf('Weighted-center angle : %8.3f deg\n', rad2deg(theta_w));
fprintf('Fit center angle      : %8.3f deg\n', rad2deg(theta_best));
fprintf('Fit sigma             : %8.3f deg\n', rad2deg(sigma_best));
fprintf('Fit FWHM              : %8.3f deg\n', rad2deg(fwhm_best));

date_center = t_ref + days(theta_best/(2*pi) * 365.25);
fprintf('Approx center date    : %s\n', datestr(date_center));





%% ESTIMATE CONE CENTER
% Choose the ESA step that best shows the cone
esa_fit = 3;

idx_fit = ram(:,2) == esa_fit;
time_fit = time(idx_fit);
rate_fit = ram(idx_fit,11);

% Sort by time
[time_fit, isrt] = sort(time_fit);
rate_fit = rate_fit(isrt);

% Orbital angle in radians, wrapped to [0, 2pi)
theta_fit_data = 2*pi * days(time_fit - t_ref) / 365.25;
theta_fit_data = mod(theta_fit_data, 2*pi);

% Remove bad values
good = isfinite(theta_fit_data) & isfinite(rate_fit);
theta_fit_data = theta_fit_data(good);
rate_fit = rate_fit(good);

% Quick weighted center estimate
% Use only the brighter part of the distribution
pct_use = 70;   % top 30% brightest points
thr = prctile(rate_fit, pct_use);
use = rate_fit >= thr;

theta_use = theta_fit_data(use);
w_use = rate_fit(use);

xw = sum(w_use .* cos(theta_use));
yw = sum(w_use .* sin(theta_use));

theta_w = atan2(yw, xw);
if theta_w < 0
    theta_w = theta_w + 2*pi;
end

% Gaussian-like fit in wrapped angle
% initial guesses
bkg0 = min(rate_fit);
[pmax, imax] = max(rate_fit);
A0 = pmax - bkg0;
theta0 = theta_w;           % use weighted center as initial guess
sigma0 = deg2rad(20);       % initial width guess

p0 = [bkg0, A0, theta0, sigma0];

opts = optimset('Display','off', 'MaxFunEvals',5000, 'MaxIter',5000);
pbest = fminsearch(@(p) cone_chisq(p, theta_fit_data, rate_fit), p0, opts);

bkg_best   = pbest(1);
A_best     = pbest(2);
theta_best = mod(pbest(3), 2*pi);
sigma_best = abs(pbest(4));
fwhm_best  = 2.355 * sigma_best;

fprintf('\n=== Cone estimate for ESA %d ===\n', esa_fit);
fprintf('Weighted-center angle : %8.3f deg\n', rad2deg(theta_w));
fprintf('Fit center angle      : %8.3f deg\n', rad2deg(theta_best));
fprintf('Fit sigma             : %8.3f deg\n', rad2deg(sigma_best));
fprintf('Fit FWHM              : %8.3f deg\n', rad2deg(fwhm_best));

date_center = t_ref + days(theta_best/(2*pi) * 365.25);
fprintf('Approx center date    : %s\n', datestr(date_center));

% Overlay gaussian FWHF on orbit plot
figure(3); hold on

r_line = 1.18;
plot([0, r_line*cos(theta_best)], [0, r_line*sin(theta_best)], ...
    'm-', 'LineWidth', 2.5, 'DisplayName',...
    sprintf('Center (ESA %d)', esa_fit));

plot(r_line*cos(theta_best), r_line*sin(theta_best), ...
    'mo', 'MarkerFaceColor','r', 'MarkerSize', 8, ...
    'HandleVisibility','off');

text(1.4*cos(theta_best), 1.4*sin(theta_best), ...
    sprintf('Center\nESA %d', esa_fit), ...
    'Color', 'm', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'FontWeight','bold');

% legend('Location','eastoutside')

% Optional: show fitted angular width on orbit plot
% draw dashed lines at center +/- sigma
theta_lo = theta_best - sigma_best;
theta_hi = theta_best + sigma_best;

plot([0, r_line*cos(theta_lo)], [0, r_line*sin(theta_lo)], ...
    'm--', 'LineWidth', 1.2, 'HandleVisibility','off');
plot([0, r_line*cos(theta_hi)], [0, r_line*sin(theta_hi)], ...
    'm--', 'LineWidth', 1.2, 'HandleVisibility','off');


%% FIT DIAGNOSTIC PLOT: rate vs orbital angle
theta_plot_deg = rad2deg(theta_fit_data);
[theta_plot_deg, isrt2] = sort(theta_plot_deg);
rate_plot = rate_fit(isrt2);

theta_model = linspace(0, 2*pi, 1000)';
rate_model = cone_model(pbest, theta_model);

figure(4); clf
plot(theta_plot_deg, rate_plot, 'ko', 'MarkerFaceColor', [0.3 0.3 0.3])
hold on
plot(rad2deg(theta_model), rate_model, 'r-', 'LineWidth', 2)

xline(rad2deg(theta_best), 'r--', 'LineWidth', 1.5, ...
    'Label', sprintf('Center = %.1f^\\circ', rad2deg(theta_best)), ...
    'LabelOrientation','horizontal');

xlabel('Orbital angle (deg)')
ylabel('Rate average')
title(sprintf('Cone fit diagnostic, ESA %d', esa_fit))
grid on
box on
xlim([0 360])

%% HELPER FUNCTIONS
function y = cone_model(p, theta)
    bkg = p(1);
    A = p(2);
    theta0 = p(3);
    sigma = abs(p(4));

    dth = angdiff_wrap(theta, theta0);
    y = bkg + A .* exp(-0.5 * (dth ./ sigma).^2);
end

function chi2 = cone_chisq(p, theta, rate)
    ymod = cone_model(p, theta);
    chi2 = nansum((rate - ymod).^2);

    % soft penalties to keep the fit sensible
    if p(2) < 0
        chi2 = chi2 + 1e8;
    end
    if abs(p(4)) > pi/2
        chi2 = chi2 + 1e8;
    end
end

function d = angdiff_wrap(a, b)
    d = atan2(sin(a-b), cos(a-b));
end



%% Function

function ang_mean = circ_mean_deg(ang_deg)
    ang_deg = ang_deg(~isnan(ang_deg));
    
    if isempty(ang_deg)
        ang_mean = NaN;
        return
    end

    ang_rad = deg2rad(ang_deg);
    x = mean(cos(ang_rad));
    y = mean(sin(ang_rad));

    ang_mean = rad2deg(atan2(y, x));
    
    if ang_mean < 0
        ang_mean = ang_mean + 360;
    end
end


function lon = get_earth_lon_eclipj2000(t_in, tls_kernel, spk_kernel)

    % Preserve original timezone if present, but convert to UTC for SPICE
    if isdatetime(t_in)
        if isempty(t_in.TimeZone)
            t_in.TimeZone = 'UTC';
        else
            t_in.TimeZone = 'UTC';
        end
    else
        error('Input must be datetime.')
    end

    % Load kernels
    cspice_kclear
    cspice_furnsh(tls_kernel)
    cspice_furnsh(spk_kernel)

    % Convert datetime -> UTC strings
    utc_str = cellstr(datestr(t_in, 'yyyy-mm-dd HH:MM:SS'));

    % UTC -> ET
    et = zeros(1, numel(utc_str));
    for ii = 1:numel(utc_str)
        et(ii) = cspice_str2et(utc_str{ii});
    end

    % Earth state relative to Sun in ECLIPJ2000
    [state, ~] = cspice_spkezr('EARTH', et, 'ECLIPJ2000', 'NONE', 'SUN');

    x = state(1,:);
    y = state(2,:);

    lon = atan2(y, x);
    lon = mod(lon, 2*pi);

    lon = lon(:);
end

function y = model_1g(p, th)
    b  = p(1);
    a  = p(2);
    mu = p(3);
    s  = abs(p(4));

    d = wrapdiff(th, mu);
    y = b + a .* exp(-0.5 * (d./s).^2);
end

function val = chisq_1g(p, th, ydat)
    ymod = model_1g(p, th);
    val = nansum((ydat - ymod).^2);

    if p(2) < 0 || abs(p(4)) > pi/2
        val = val + 1e8;
    end
end

function y = model_2g(p, th)
    b   = p(1);
    a1  = p(2);
    mu1 = p(3);
    s1  = abs(p(4));
    a2  = p(5);
    mu2 = p(6);
    s2  = abs(p(7));

    d1 = wrapdiff(th, mu1);
    d2 = wrapdiff(th, mu2);

    y = b ...
      + a1 .* exp(-0.5 * (d1./s1).^2) ...
      + a2 .* exp(-0.5 * (d2./s2).^2);
end

function val = chisq_2g(p, th, ydat)
    ymod = model_2g(p, th);
    val = nansum((ydat - ymod).^2);

    if p(2) < 0 || p(5) < 0 || abs(p(4)) > pi/2 || abs(p(7)) > pi/2
        val = val + 1e8;
    end
end

function d = wrapdiff(a, b)
    d = atan2(sin(a-b), cos(a-b));
end

function y = model_3g(p, th)
    b   = p(1);

    a1  = p(2);  mu1 = p(3);  s1 = abs(p(4));
    a2  = p(5);  mu2 = p(6);  s2 = abs(p(7));
    a3  = p(8);  mu3 = p(9);  s3 = abs(p(10));

    d1 = wrapdiff(th, mu1);
    d2 = wrapdiff(th, mu2);
    d3 = wrapdiff(th, mu3);

    y = b ...
      + a1 .* exp(-0.5 * (d1./s1).^2) ...
      + a2 .* exp(-0.5 * (d2./s2).^2) ...
      + a3 .* exp(-0.5 * (d3./s3).^2);
end

function val = chisq_3g(p, th, ydat)
    ymod = model_3g(p, th);
    val = nansum((ydat - ymod).^2);

    % simple penalties
    if p(2) < 0 || p(5) < 0 || p(8) < 0
        val = val + 1e8;
    end
    if abs(p(4)) > pi/2 || abs(p(7)) > pi/2 || abs(p(10)) > pi/2
        val = val + 1e8;
    end
end

% function d = wrapdiff(a, b)
%     d = atan2(sin(a-b), cos(a-b));
% end