%% MAIN

close all; clearvars; clc; format long; warning('off','all')

fld = '/Users/mitchellshen/CAVA_IDL/IMAP-Lo June 10th 2026 Maps/';
qmp = 'l3/2025/11/';
% fld = '/Users/mitchellshen/Library/CloudStorage/Dropbox-Personal/';
% qmp = 'IMAP_ILO/SCI_ENA/maps_';
% qmp = 'IMAP_ILO/SCI_ENA/p90maps_arxiv/';
% qmp = 'IMAP_ILO/SCI_ENA/PAC_correction/hydrogen_map/outdir/pivot_90/maps/';

%==================== USER SWITCH ====================%
flip_view = true;   % true = outside-looking-in style
slice_view= false;
boot      = false;
sput      = false;
rate      = false;
save      = false;
species   = 'H';
pivot     = 90;
esa_step  = 3;
esa_volt  = {'16','30','56','106','200','404','787'};

%==================== READ DATA ====================%
E = [];
if boot && sput
    error('boot and sput cannot both be true at the same time.');
end
if boot
    % e.g., map_flux_7_Hy_boot_cor
    qmp = 'IMAP_ILO/SCI_ENA/p90maps_BScorr/';
    Z = table2array(readtable([fld qmp 'map_flux_' num2str(esa_step) '_Hy_boot_cor.csv']));
elseif sput
    qmp = 'IMAP_ILO/SCI_ENA/p90maps_BScorr/';
    Z = table2array(readtable([fld qmp 'map_flux_' num2str(esa_step) '_Hy_sput_cor.csv']));
elseif rate
     % e.g., map_flux_esa7
    qmp = [qmp species '_p' num2str(pivot,'%03d') '/'];
    Z = table2array(readtable([fld qmp 'map_rate_esa' num2str(esa_step) '.csv']));
else
    % e.g., map_flux_esa7
    qmp = [qmp species '_p' num2str(pivot,'%03d') '/'];
    Z = table2array(readtable([fld qmp 'map_flux_esa' num2str(esa_step) '.csv']));
    E = table2array(readtable([fld qmp 'map_func_esa' num2str(esa_step) '.csv']));
end
Z(1,:) = [];
if ~isempty(E)
    E(1,:) = [];
end

lat =  -90:6: 90-6;
lon = (0+6:6:360);

%==================== LONGITUDE SETUP ====================%
% Convert 0:354 --> [-180,180)
lon = mod(lon + 180, 360) - 180;

lat_center = 5;
lon_center = 105;

% Flip data/view if desired
if flip_view
    lon = -lon;
end

% Sort longitude and reorder data
[lon, isrt] = sort(lon);
Z = Z(:, isrt);
if ~isempty(E)
    E = E(:, isrt);
end

%==================== MASK NO-DATA ====================%
Zmask = ~isnan(Z) & (Z > 0);
Z(~Zmask) = NaN;
if ~isempty(E)
    E(~Zmask) = NaN;
end

[Lon, Lat] = meshgrid(lon, lat);

%==================== FIGURE ====================%
figure(1);
set(gcf, 'Color', [1 1 1]);
ax = axesm('mollweid', ...
    'Origin',[lat_center lon_center 0], ...
    'MapLatLimit', [-90 90], ...
    'MapLonLimit', [-180 180], ...
    'Frame', 'on', ...
    'Grid', 'on', ...
    'MeridianLabel', 'off', ...   % turn OFF automatic labels
    'ParallelLabel', 'off', ...
    'MLineLocation', 30, ...
    'PLineLocation', 30, ...
    'GLineWidth', 1.0, ...
    'Gcolor', [1 0.6 0], ...
    'FEdgeColor', [0 0 0.5], ...
    'FontSize', 12);

% mlabel on
% setm(ax, 'MLabelParallel', 0, 'FontColor', [1 1 1])

axis off
tightmap

% White outside globe, dark inside globe
set(ax, 'Color', 'none')
% setm(ax, 'FFaceColor', [0.25 0 0.35]) % dark purple
setm(ax, 'FFaceColor', [0 0 0])


%==================== PLOT MAP ====================%
% Wrap-pad longitude (add lon(1)+360 as an extra column, duplicating first column's data)
lon_ext = [lon, lon(1) + 360];
Z_ext   = [Z,   Z(:,1)];
if ~isempty(E)
    E_ext = [E, E(:,1)];
end

% Pad latitude similarly (duplicate last row so the top band isn't dropped)
lat_ext = [lat, lat(end) + (lat(2)-lat(1))];
Z_ext   = [Z_ext; Z_ext(end,:)];
if ~isempty(E)
    E_ext = [E_ext; E_ext(end,:)];
end

[Lon, Lat] = meshgrid(lon_ext, lat_ext);
Zmask = ~isnan(Z_ext) & (Z_ext > 0);
Z_ext(~Zmask) = NaN;

% h = surfm(Lat, Lon, Z);
h = surfm(Lat, Lon, Z_ext);
set(h, ...
    'EdgeColor', 'none', ...
    'FaceColor', 'flat', ...
    'FaceAlpha', 'flat', ...
    'AlphaData', double(Zmask), ...
    'AlphaDataMapping', 'none');

set(gca, 'ColorScale', 'log')
colormap("turbo")

title(['Pivot:' num2str(pivot)], 'FontSize', 18);

%==================== MANUAL MERIDIAN LABELS ====================%
lablat = 0;   % latitude where labels sit (equator)

meridians = -180:30:180;

for k = 1:numel(meridians)
    m = meridians(k);

    % Flip label meaning relative to plotted meridian
    val = -m;

    % Wrap into [-180,180]
    val = mod(val + 180, 360) - 180;

    % Format as signed degree (recommended)
    if abs(val) == 180
        txt = '180°';
    elseif val == 0
        txt = '0°';
    else
        txt = sprintf('%d°', val);
    end

    textm(lablat, m, txt, ...
        'Color', [1 1 1], ...
        'FontSize', 12, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle');
end

%==================== COLORBAR ====================%
cb = colorbar('southoutside');
if ~rate
cb.Label.String = ['Intensity (counts cm^{-2} s^{-1} sr^{-1} keV^{-1})' ...
    ' [ECLIPJ2000] at ESA (eV) = ' char(esa_volt(esa_step))];
else
cb.Label.String = ['Rate (counts s^{-1})' ...
    ' [ECLIPJ2000] at ESA (eV) = ' char(esa_volt(esa_step))];
end
cb.FontSize = 14;

if ~rate
switch esa_step
    case 1
        clim([1e3 1e6])
    case {2,3,4}
        clim([1e2 1e6])
    case 5
        clim([1e1 1e4])
    case 6
        clim([0 500])
        % set(gca, 'ColorScale', 'lin')
    case 7
        clim([1e1  1000])
end
else
    clim([1e-3  1e1])
end


%==================== IBEX RIBBON GUIDE ====================%

% Funsten et al. 2010 centered on (221, 39)
% and is 36 degrees wide (radius 71.6 degrees with +/- 18 degrees around it). 
% It is for illustrations to guide the eye on ribbon estimate boundaries only
% We also used it in McComas papers for 5-yr and 7-yr of IBEX observations.

% Later on, an updated trace used the weighted center from 
% different Ribbon fits instead, averaged over all energies (0.7–4.3 keV) 
% and all nine years of IBEX-Hi data (2009–2017), 
% from our 2019 paper on the ribbon variations. These are:
% Center ( b_lon_ecl = 218.33, b_lat_ecl = 40.38) and radius=74.81



% Ribbon center (Funsten 2010)
lon0  = 218.33; % 221.00;
lat0  =  40.38; %  39.00;
width = 36;

% ribbon radii
r_mid =  74.81; %  71.6;
r_in  = r_mid - width./2;
r_out = r_mid + width./2;

% angle around circle
az = linspace(0,360,721);

% generate circle points
[lat_mid, lon_mid] = smallcircle_deg(lat0, lon0, r_mid, az);
[lat_in , lon_in ] = smallcircle_deg(lat0, lon0, r_in , az);
[lat_out, lon_out] = smallcircle_deg(lat0, lon0, r_out, az);

% Voyager dots
lon_v1 = mod(255.0 + 180, 360) - 180;
lat_v1 =  35.0;
lon_v2 = mod(290.4 + 180, 360) - 180;
lat_v2 = -32.2;

% Nose and Tail
lon_nose = mod(259.200622 + 180, 360) - 180;
lat_nose =  5.116296;
lon_tail = mod( 79.200622 + 180, 360) - 180;
lat_tail = -5.116296;

% lons = np.array([259.200622, 255.499300, 289.800000, 79.200622, 221,150.,0.0, 41, 250.624523, 70.624523])
% lats = np.array([5.116296, 34.999941, -35.600000, -5.116296, 39, 0.0, 0.0, -39, -41.242979, 41.242979])
% labels=[ 'Nose', 'V1','V2','Tail', 'RibC','Stbd', 'Port','DwnB','ULSR', 'DLSR']


% 🔥 IMPORTANT: match your map convention
if flip_view
    lon_mid = -lon_mid;
    lon_in  = -lon_in;
    lon_out = -lon_out;
    lon_v1  = -lon_v1;
    lon_v2  = -lon_v2;
    lon_nose= -lon_nose;
    lon_tail= -lon_tail;
    lon0    = -lon0;
end

% wrap again after flipping
lon_mid = mod(lon_mid + 180, 360) - 180;
lon_in  = mod(lon_in  + 180, 360) - 180;
lon_out = mod(lon_out + 180, 360) - 180;
lon_v1  = mod(lon_v1  + 180, 360) - 180;
lon_v2  = mod(lon_v2  + 180, 360) - 180;
lon_nose= mod(lon_nose  + 180, 360) - 180;
lon_tail= mod(lon_tail  + 180, 360) - 180;
lon0    = mod(lon0    + 180, 360) - 180;

% plot ribbon
plotm(lat_mid, lon_mid, 'w-',  'LineWidth', 2)
plotm(lat_in , lon_in , 'w--', 'LineWidth', 1.5)
plotm(lat_out, lon_out, 'w--', 'LineWidth', 1.5)


% [lat_mid, lon_mid] = smallcircle_deg(lat0, lon0, r_mid, az); NO NEED

%==================== RIBBON  MARKERS ====================%

plotm(lat0, lon0, 'wo', ...
    'MarkerFaceColor','w', 'MarkerSize',6)
textm(lat0, lon0, 'RC', ...
    'Color','w','FontSize',12,'FontWeight','bold', ...
    'HorizontalAlignment','right','VerticalAlignment','bottom')

%==================== VOYAGER MARKERS ====================%

plotm(lat_v1, lon_v1, 'wo', ...
    'MarkerFaceColor','w', 'MarkerSize',6)
textm(lat_v1, lon_v1, '  V1', ...
    'Color','w','FontSize',12,'FontWeight','bold',...
    'HorizontalAlignment','right', 'VerticalAlignment','bottom')

plotm(lat_v2, lon_v2, 'wo', ...
    'MarkerFaceColor','w', 'MarkerSize',6)
textm(lat_v2, lon_v2, '  V2', ...
    'Color','w','FontSize',12,'FontWeight','bold',...
    'HorizontalAlignment','right', 'VerticalAlignment','top')

%==================== Nose/Tail MARKERS ====================%
plotm(lat_nose, lon_nose, 'wo', ...
    'MarkerFaceColor','w', 'MarkerSize',6)
textm(lat_nose, lon_nose-3, '  Nose', ...
    'Color','w','FontSize',12,'FontWeight','bold',...
    'HorizontalAlignment','left', 'VerticalAlignment','bottom')
   %'HorizontalAlignment','right', 'VerticalAlignment','bottom')

plotm(lat_tail, lon_tail, 'wo', ...
    'MarkerFaceColor','w', 'MarkerSize',6)
textm(lat_tail, lon_tail, '  Tail', ...
    'Color','w','FontSize',12,'FontWeight','bold',...
    'HorizontalAlignment','right', 'VerticalAlignment','top')

%==================== CROSS-RIBBON CUTS ====================%

% keep everything except lon in [-140, -90]
use = ~(lon_mid >= 60 & lon_mid <= 130) & ...   
      (lat_mid >= -90 & lat_mid <= 90);
      % ~(lon_mid >= -180 & lon_mid <= -140) &


idx_full = find(use);

% keep only points that have neighbors for tangent estimation
idx_full = idx_full(idx_full > 1 & idx_full < numel(lon_mid));

% optional downsampling after selection
s_cnt = 16;
idx_use = idx_full(1:s_cnt:end);

% cut distances in degrees across ribbon
% s_cut = lat0-90-12-3:6:lat0+90-24-3;   % inward/outward from ribbon centerline
s_cut = lat0-75:6:lat0; 

all_cut_lat = cell(size(idx_use));
all_cut_lon = cell(size(idx_use));

for ii = 1:numel(idx_use)
    i = idx_use(ii);

    % neighboring points to estimate tangent
    p_prev = sph2cart_unit(lon_mid(i-1), lat_mid(i-1));
    p_here = sph2cart_unit(lon_mid(i),   lat_mid(i));
    p_next = sph2cart_unit(lon_mid(i+1), lat_mid(i+1));

    % local tangent along ribbon
    t = p_next - p_prev;
    t = t / norm(t);

    % local radial direction on sphere
    r = p_here / norm(p_here);

    % local normal within tangent plane, perpendicular to ribbon
    n = cross(r, t);
    n = n / norm(n);

    % build cut by rotating p_here around axis = t
    lat_cut = nan(size(s_cut));
    lon_cut = nan(size(s_cut));

    for j = 1:numel(s_cut)
        p_cut = rotate_about_axis(p_here, t, deg2rad(s_cut(j)));
        [lon_cut(j), lat_cut(j)] = cart2sph_deg(p_cut);
    end

    all_cut_lat{ii} = lat_cut;
    all_cut_lon{ii} = lon_cut;

    % overlay cut on map
    if slice_view
        plotm(lat_cut, lon_cut, 'm-', 'LineWidth', 1.5)
    end
end

%

%==================== SAMPLE DATA ALONG CUTS ====================%
Fz = griddedInterpolant({lat, lon}, Z, 'linear', 'none');

if ~isempty(E)
    Fe = griddedInterpolant({lat, lon}, E, 'linear', 'none');
end

all_cut_vals = cell(size(idx_use));
all_cut_errs = cell(size(idx_use));

for ii = 1:numel(idx_use)
    lat_cut = all_cut_lat{ii};
    lon_cut = all_cut_lon{ii};

    vals_cut = Fz(lat_cut, lon_cut);
    all_cut_vals{ii} = vals_cut;

    if ~isempty(E)
        errs_cut = Fe(lat_cut, lon_cut);
        all_cut_errs{ii} = errs_cut;
    end
end

figure(2)
hold on
for ii = 1:numel(all_cut_vals)
    plot(s_cut, all_cut_vals{ii}, '-','LineWidth',2)
end
grid on; box on;
set(gca,'YScale','log')
xlabel('Angular distance across ribbon (deg)')
ylabel('Intensity')

%
%==================== SUM PROFILE ====================%

P = nan(numel(all_cut_vals), numel(s_cut));
Pe = nan(numel(all_cut_errs), numel(s_cut));

for ii = 1:numel(all_cut_vals)
    P(ii,:) = all_cut_vals{ii};
    if ~isempty(E)
        Pe(ii,:) = all_cut_errs{ii};
    end
end

sum_profile =  nansum(P, 1);   % sum across cuts
ave_profile = nanmean(P, 1);   % ave across cuts
ave_profile_err = sqrt(nansum(Pe.^2, 1)) ./ sum(~isnan(P), 1);
sum_profile_err = sqrt(nansum(Pe.^2, 1));

figure(3)
plot(s_cut, ave_profile, 'k-', 'LineWidth', 2)
set(gca,'YScale','log')
xlabel('Angular distance across ribbon (deg)')
ylabel('Average Intensity')
title('Stacked Ribbon Profile (Average)')
grid on

figure(4)
errorbar(s_cut, ave_profile, ave_profile_err, 'k-', 'LineWidth', 1.5)
set(gca,'YScale','log')
xlabel('Angular distance across ribbon (deg)')
ylabel('Average Intensity')
title('Stacked Ribbon Profile (Average ± 1\sigma)')
grid on

%
P = nan(numel(idx_use), numel(s_cut));
for ii = 1:numel(idx_use)
    P(ii,:) = all_cut_vals{ii};
end

figure(5)
imagesc(1:numel(idx_use), s_cut, P')
set(gca,'YDir','reverse')
colormap("turbo")
colorbar
xlabel('Cut index along ribbon')
ylabel('Angular distance across ribbon (deg)')
set(gca, 'ColorScale', 'log')

figure(1)

if save
    T = table(s_cut(:), ave_profile(:), ave_profile_err(:),...
        'VariableNames', {'s_cut_deg','average_intensity','ave_profile_err'});
    writetable(T, [fld qmp 'map_ave_flux_esa' num2str(esa_step) '.csv']);
else
end

%% Profiles
% 1. File path
file = '/Users/mitchellshen/Library/CloudStorage/Dropbox-Personal/IMAP_ILO/SCI_ENA/p90maps_ribbon/map_ave_flux_combo.xlsx';

% 2. Import as table (robust)
opts = detectImportOptions(file);
opts.VariableNamingRule = 'preserve';   % keep original column names
T = readtable(file, opts);

% 3. Inspect data
disp('Column names:')
disp(T.Properties.VariableNames)

disp('First few rows:')
head(T)

% 4. Extract columns (EDIT names if needed)
% >>> CHANGE these names based on your actual column headers <<<
% s_cut       = T.s_cut;         % e.g., angular distance
% ave_profile = T.ave_profile;   % e.g., averaged intensity

% Ensure column vectors
M = table2array(T);

figure('Color',[1 1 1])   % FIXED version of your earlier error
for i = 2: 8 %size(M,2)
    plot(M(:,1), M(:,i), 'LineWidth', 2); hold on;
end
% for i = 9: 15 %size(M,2)
%     plot(M(:,1), M(:,i), 'LineStyle', '-.', 'LineWidth', 2); hold on;
% end
set(gca,'YScale','log')
xlabel('Angular distance across ribbon (deg)')
ylabel('Average Intensity')
title('Stacked Ribbon Profile (Average)')
grid on
box on

% 6. Export to CSV
outTable = table(s_cut, ave_profile);
writetable(outTable, 'ribbon_profile.csv');

disp('CSV exported: ribbon_profile.csv')

% 7. OPTIONAL: If file is actually a 2D flux map
% Uncomment if your Excel is a grid instead of columns
% M = readmatrix(file);
% figure
% imagesc(M)
% set(gca,'YDir','normal')
% colorbar
% title('ENA Flux Map')

%% WOWOW
x = x_all;
n_esa = size(Y_resid, 2);

A_rib_fit   = nan(1, n_esa);
mu_rib_fit  = nan(1, n_esa);
sig_rib_fit = nan(1, n_esa);

Y_ribbon_fit = nan(size(Y_resid));

mu_fix = 0;
fit_half_width = 22;

sigma_grid = linspace(6, 15, 200);

for j = 1:n_esa

    yres = Y_resid(:,j);

    good = isfinite(x) & isfinite(yres) & (yres > 0);
    xg = x(good);
    yg = yres(good);

    fit_mask = xg >= -fit_half_width & xg <= fit_half_width;
    x_fit = xg(fit_mask);
    y_fit = yg(fit_mask);

    if numel(x_fit) < 5
        warning('ESA%d residual does not have enough points for ribbon fitting.', j);
        continue
    end

    best_sse = inf;
    best_A = nan;
    best_sigma = nan;

    for k = 1:numel(sigma_grid)
        sigma_try = sigma_grid(k);

        g = exp(-((x_fit - mu_fix).^2) ./ (2*sigma_try^2));

        A_try = sum(y_fit .* g) / sum(g.^2);

        y_try = A_try .* g;

        sse = sum((y_fit - y_try).^2);

        if sse < best_sse
            best_sse = sse;
            best_A = A_try;
            best_sigma = sigma_try;
        end
    end

    gfull = exp(-((x - mu_fix).^2) ./ (2*best_sigma^2));

    y_rib_full = nan(size(x));
    local_plot_mask = abs(x - mu_fix) <= 25;
    y_rib_full(local_plot_mask) = best_A .* gfull(local_plot_mask);

    A_rib_fit(j)   = best_A;
    mu_rib_fit(j)  = mu_fix;
    sig_rib_fit(j) = best_sigma;

    Y_ribbon_fit(:,j) = y_rib_full;
end

disp('Fixed-center, bounded-width ribbon proxy Gaussian parameters:')
for j = 1:n_esa
    fprintf('ESA%d: A = %.4g, mu = %.3f deg, sigma = %.3f deg\n', ...
        j, A_rib_fit(j), mu_rib_fit(j), sig_rib_fit(j));
end

clr = lines(n_esa);

figure('Color','w')
for j = 1:n_esa
    yres = Y_resid(:,j);
    yrib = Y_ribbon_fit(:,j);

    good_res = isfinite(x) & isfinite(yres) & yres > 0;
    good_fit = isfinite(x) & isfinite(yrib) & yrib > 0;

    semilogy(x(good_res), yres(good_res), '.', 'Color', clr(j,:), 'MarkerSize', 14); hold on
    semilogy(x(good_fit), yrib(good_fit), '-', 'Color', clr(j,:), 'LineWidth', 2.5);
end
xlabel('Slice deg')
ylabel('Residual after ISN subtraction')
title('Fixed-center, bounded-width ribbon proxy Gaussian fits')
grid on
box on
xlim([-90 110])
legend('ESA1 resid','ESA1 rib fit','ESA2 resid','ESA2 rib fit','ESA3 resid','ESA3 rib fit', ...
       'ESA4 resid','ESA4 rib fit','ESA5 resid','ESA5 rib fit','ESA6 resid','ESA6 rib fit', ...
       'ESA7 resid','ESA7 rib fit','Location','eastoutside')


%% WOW

x_all = M(:,1);
Y_all = M(:,2:8);

n_esa = size(Y_all, 2);

A_fit     = nan(1, n_esa);
mu_fit    = nan(1, n_esa);
sigma_fit = nan(1, n_esa);

Y_isn   = nan(size(Y_all));
Y_resid = nan(size(Y_all));

gauss_fun = @(p, xx) p(1) .* exp(-((xx - p(2)).^2) ./ (2*p(3)^2));

for j = 1:n_esa

    x = x_all;
    y = Y_all(:,j);

    good = isfinite(x) & isfinite(y) & y > 0;
    xg = x(good);
    yg = y(good);

    fit_mask = xg >= -65 & xg <= -10;
    x_fit = xg(fit_mask);
    y_fit = yg(fit_mask);

    if numel(x_fit) < 5
        warning('ESA%d does not have enough points for fitting.', j);
        continue
    end

    A0 = max(y_fit);
    mu0 = -33;
    sigma0 = 12;
    p0 = [A0, mu0, sigma0];

    objfun = @(p) ...
        sum((y_fit - gauss_fun([p(1), p(2), abs(p(3))], x_fit)).^2) ...
        + 1e6 * (p(2) + 33).^2;

    p_best = fminsearch(objfun, p0);

    A     = p_best(1);
    mu    = p_best(2);
    sigma = abs(p_best(3));

    y_isn_full = nan(size(y));
    y_isn_full(good) = gauss_fun([A, mu, sigma], xg);

    y_resid_full = y - y_isn_full;

    A_fit(j) = A;
    mu_fit(j) = mu;
    sigma_fit(j) = sigma;

    Y_isn(:,j) = y_isn_full;
    Y_resid(:,j) = y_resid_full;
end

disp('Fitted Gaussian parameters:')
for j = 1:n_esa
    fprintf('ESA%d: A = %.4g, mu = %.3f deg, sigma = %.3f deg\n', ...
        j, A_fit(j), mu_fit(j), sigma_fit(j));
end

clr = lines(n_esa);

figure('Color','w')
for j = 1:n_esa
    y = Y_all(:,j);
    yfit = Y_isn(:,j);

    good_data = isfinite(x_all) & isfinite(y) & y > 0;
    good_fit  = isfinite(x_all) & isfinite(yfit) & yfit > 0;

    semilogy(x_all(good_data), y(good_data), '.', 'Color', clr(j,:), 'MarkerSize', 14); hold on
    semilogy(x_all(good_fit),  yfit(good_fit), '-', 'Color', clr(j,:), 'LineWidth', 2.5);
end
xlabel('Slice deg')
ylabel('Signal')
title('All ESA curves with ISN Gaussian fits')
grid on
box on
xlim([-90 110])
legend('ESA1 data','ESA1 fit','ESA2 data','ESA2 fit','ESA3 data','ESA3 fit', ...
       'ESA4 data','ESA4 fit','ESA5 data','ESA5 fit','ESA6 data','ESA6 fit', ...
       'ESA7 data','ESA7 fit','Location','eastoutside')

figure('Color','w')
for j = 1:n_esa
    yres = Y_resid(:,j);
    good = isfinite(x_all) & isfinite(yres);
    plot(x_all(good), yres(good), '-', 'Color', clr(j,:), 'LineWidth', 1.8); hold on
end
yline(0,'k--','LineWidth',1.2)
xlabel('Slice deg')
ylabel('Residual = Data - ISN Gaussian')
title('Residuals for all ESA curves')
grid on
box on
xlim([-90 110])
legend('ESA1','ESA2','ESA3','ESA4','ESA5','ESA6','ESA7','Location','eastoutside')

figure('Color','w')
for j = 1:n_esa
    yres = Y_resid(:,j);
    good = isfinite(x_all) & isfinite(yres) & yres > 0;
    semilogy(x_all(good), yres(good), '-', 'Color', clr(j,:), 'LineWidth', 2); hold on
end
xlabel('Slice deg')
ylabel('Positive residual')
title('Positive residuals after ISN subtraction')
grid on
box on
xlim([-90 110])
legend('ESA1','ESA2','ESA3','ESA4','ESA5','ESA6','ESA7','Location','eastoutside')

%%
% Assume:
% M(:,1) = slice angle (deg)
% M(:,5) = ESA4

x = M(:,1);
y = M(:,5);

% Keep valid positive data
good = isfinite(x) & isfinite(y) & y > 0;
x = x(good);
y = y(good);

% Select fitting region (ISN-dominated side)
fit_mask = x >= -65 & x <= -10;
x_fit = x(fit_mask);
y_fit = y(fit_mask);

% Gaussian model: y = A * exp(-(x-mu)^2/(2*sigma^2))
gauss_fun = @(p, xx) p(1) .* exp(-((xx - p(2)).^2) ./ (2*p(3)^2));

% Initial guess (center near -33 deg)
A0     = max(y_fit);
mu0    = -33;
sigma0 = 12;
p0 = [A0, mu0, sigma0];

% Objective function with penalty to keep center near -33
objfun = @(p) ...
    sum((y_fit - gauss_fun([p(1), p(2), abs(p(3))], x_fit)).^2) ...
    + 1e6 * (p(2) + 33).^2;

% Fit
p_best = fminsearch(objfun, p0);

A     = p_best(1);
mu    = p_best(2);
sigma = abs(p_best(3));

% Evaluate Gaussian over full range
y_isn = gauss_fun([A, mu, sigma], x);

% Residual (subtract ISN component)
y_resid = y - y_isn;

% Plot ESA4 and Gaussian fit
figure('Color','w')
semilogy(x, y, 'kx-', 'LineWidth', 1.5, 'MarkerSize', 8); hold on
semilogy(x, y_isn, 'm-', 'LineWidth', 3)
xlabel('Slice deg')
ylabel('ESA4')
title('ESA4 with ISN Gaussian Fit')
legend('ESA4 data', 'ISN Gaussian fit', 'Location', 'best')
grid on
box on
xlim([-90 110])

% Plot residual (linear scale)
figure('Color','w')
plot(x, y_resid, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6); hold on
yline(0, 'k--', 'LineWidth', 1.2)
xlabel('Slice deg')
ylabel('Residual = ESA4 - ISN Gaussian')
title('Residual after subtracting ISN contribution')
grid on
box on
xlim([-90 110])

% Plot everything together (log scale, positive residual only)
figure('Color','w')
semilogy(x, y, 'k.-', 'LineWidth', 1.2, 'MarkerSize', 14); hold on
semilogy(x, y_isn, 'm-', 'LineWidth', 3)

mask_pos = y_resid > 0;
semilogy(x(mask_pos), y_resid(mask_pos), 'r.-', 'LineWidth', 1.5, 'MarkerSize', 14)

xlabel('Slice deg')
ylabel('Signal')
title('ESA4, ISN fit, and residual (positive only)')
legend('ESA4', 'ISN Gaussian', 'Residual (>0)', 'Location', 'best')
grid on
box on
xlim([-90 110])

% Print fitted parameters
fprintf('ISN Gaussian fit parameters for ESA4:\n');
fprintf('A     = %.4g\n', A);
fprintf('mu    = %.3f deg\n', mu);
fprintf('sigma = %.3f deg\n', sigma);

%% HRC
%==================== RIBBON-CENTERED POLAR SKY MAP ====================%

% Ribbon center (Funsten et al. 2010)
lon0 = 221;
lat0 = 39;

% wrap to [-180,180]
lon0 = mod(lon0 + 180, 360) - 180;

% match your map convention
if flip_view
    lon0 = -lon0;
end
lon0 = mod(lon0 + 180, 360) - 180;

%--- flatten map ---%
lonv = Lon(:);
latv = Lat(:);
Zv   = Z(:);

good = ~isnan(Zv) & Zv > 0;
lonv = lonv(good);
latv = latv(good);
Zv   = Zv(good);

%==================== CORE TRANSFORM ====================%

% radial distance
rho = acosd( sind(latv).*sind(lat0) + ...
             cosd(latv).*cosd(lat0).*cosd(lonv - lon0) );

% along-ribbon angle
phi = lonv - lon0;
phi = mod(phi + 180, 360) - 180;
phi = -phi;


%==================== BIN TO GRID ====================%

phi_edges = -180:6:180;
rho_edges = 0:6:180;

phi_cent = (phi_edges(1:end-1) + phi_edges(2:end))/2;
rho_cent = (rho_edges(1:end-1) + rho_edges(2:end))/2;

Zpol = nan(numel(rho_cent), numel(phi_cent));

ibin_phi = discretize(phi, phi_edges);
ibin_rho = discretize(rho, rho_edges);

for ir = 1:numel(rho_cent)
    for ip = 1:numel(phi_cent)
        use = (ibin_rho == ir) & (ibin_phi == ip);
        if any(use)
            Zpol(ir, ip) = median(Zv(use), 'omitnan');
        end
    end
end

%==================== PLOT ====================%

[PhiGrid, RhoGrid] = meshgrid(phi_cent, rho_cent);

X = RhoGrid .* cosd(PhiGrid);
Y = RhoGrid .* sind(PhiGrid);

figure
pcolor(X, Y, Zpol)
shading flat
set(gca,'ColorScale','log')
colormap(turbo)

axis equal
axis([-185 185 -185 185])
axis off
hold on

% outer boundary
th = linspace(0,360,721);
plot(180*cosd(th), 180*sind(th), 'k-', 'LineWidth', 1.5)

% radial circles
for rr = 45:45:135
    plot(rr*cosd(th), rr*sind(th), 'k--', 'LineWidth', 0.8)
end

% spokes
for aa = 0:30:330
    plot([0 180*cosd(aa)], [0 180*sind(aa)], 'k--', 'LineWidth', 0.8)
end

% azimuth labels
for aa = 0:30:330
    xt = 200*cosd(aa);
    yt = 200*sind(aa);
    text(xt, yt, sprintf('%d', aa), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'FontSize',16);
end

% ribbon guide
r_mid = 71.6;
plot(r_mid*cosd(th), r_mid*sind(th), 'k-', 'LineWidth', 2)

%==================== V1 / V2 / RC ====================%

% Voyager directions
lon_v1 = 255.0; lat_v1 = 35.0;
lon_v2 = 290.4; lat_v2 = -32.2;

% wrap
lon_v1 = mod(lon_v1 + 180, 360) - 180;
lon_v2 = mod(lon_v2 + 180, 360) - 180;

% match map convention
if flip_view
    lon_v1 = -lon_v1;
    lon_v2 = -lon_v2;
end
lon_v1 = mod(lon_v1 + 180, 360) - 180;
lon_v2 = mod(lon_v2 + 180, 360) - 180;

% --- RC ---
plot(0, 0, 'mo', 'MarkerFaceColor','w', 'MarkerSize',8)
text(0, 0, ' RC', ...
    'Color','m','FontSize',14,'FontWeight','bold', ...
    'HorizontalAlignment','left','VerticalAlignment','bottom')

% --- V1 ---
rho_v1 = acosd( sind(lat_v1).*sind(lat0) + ...
                cosd(lat_v1).*cosd(lat0).*cosd(lon_v1 - lon0) );

phi_v1 = lon_v1 - lon0;
phi_v1 = mod(phi_v1 + 180, 360) - 180;
phi_v1 = -phi_v1;

x_v1 = rho_v1 * cosd(phi_v1);
y_v1 = rho_v1 * sind(phi_v1);

plot(x_v1, y_v1, 'mo', 'MarkerFaceColor','w', 'MarkerSize',8)
text(x_v1, y_v1, ' V1', ...
    'Color','m','FontSize',14,'FontWeight','bold', ...
    'HorizontalAlignment','left','VerticalAlignment','bottom')

% --- V2 ---
rho_v2 = acosd( sind(lat_v2).*sind(lat0) + ...
                cosd(lat_v2).*cosd(lat0).*cosd(lon_v2 - lon0) );

phi_v2 = lon_v2 - lon0;
phi_v2 = mod(phi_v2 + 180, 360) - 180;
phi_v2 = -phi_v2;

x_v2 = rho_v2 * cosd(phi_v2);
y_v2 = rho_v2 * sind(phi_v2);

plot(x_v2, y_v2, 'mo', 'MarkerFaceColor','w', 'MarkerSize',8)
text(x_v2, y_v2, ' V2', ...
    'Color','m','FontSize',14,'FontWeight','bold', ...
    'HorizontalAlignment','left','VerticalAlignment','bottom')

%==================== COLORBAR ====================%

title('Ribbon-centered polar sky map', 'FontSize', 18)

cb = colorbar;
cb.Label.Interpreter = 'none';
cb.Label.String = sprintf( ...
    'Intensity (counts cm^{-2} s^{-1} sr^{-1} keV^{-1}) [ECLIPJ2000] at ESA (eV) = %.0f', ...
    char(esa_volt(esa_step)));
cb.FontSize = 14;




%% Functions

function [lat2, lon2] = smallcircle_deg(lat1, lon1, radius_deg, az_deg)
% Generate a small circle on a sphere (degrees)

lat2 = asind( sind(lat1).*cosd(radius_deg) + ...
              cosd(lat1).*sind(radius_deg).*cosd(az_deg) );

lon2 = lon1 + atan2d( sind(az_deg).*sind(radius_deg).*cosd(lat1), ...
                      cosd(radius_deg) - sind(lat1).*sind(lat2) );

% wrap to [-180,180]
lon2 = mod(lon2 + 180, 360) - 180;
end

%
function v = sph2cart_unit(lon_deg, lat_deg)
lon = deg2rad(lon_deg);
lat = deg2rad(lat_deg);
v = [cos(lat).*cos(lon), cos(lat).*sin(lon), sin(lat)];
end

function [lon_deg, lat_deg] = cart2sph_deg(v)
v = v ./ norm(v);
lat_deg = asind(v(3));
lon_deg = atan2d(v(2), v(1));
lon_deg = mod(lon_deg + 180, 360) - 180;
end

function v_rot = rotate_about_axis(v, k, theta)
% Rodrigues rotation
k = k / norm(k);
v_rot = v*cos(theta) + cross(k, v)*sin(theta) + k*dot(k, v)*(1-cos(theta));
v_rot = v_rot / norm(v_rot);
end

function d = angdist_deg(lat1, lon1, lat2, lon2)
d = acosd( sind(lat1).*sind(lat2) + cosd(lat1).*cosd(lat2).*cosd(lon1-lon2) );
end

function phi = azimuth_about_center(lat, lon, lat0, lon0)
% azimuth around the center, measured in the local tangent plane

rc = sph2cart_unitvec(lon0, lat0);

east  = [-sind(lon0), cosd(lon0), 0];
north = [-sind(lat0).*cosd(lon0), -sind(lat0).*sind(lon0), cosd(lat0)];

v = sph2cart_unitvec(lon, lat);
proj = v - (v*rc.') .* rc;

x = proj*east.';
y = proj*north.';

phi = atan2d(y, x);
end

function v = sph2cart_unitvec(lon_deg, lat_deg)
lon = deg2rad(lon_deg);
lat = deg2rad(lat_deg);

x = cos(lat).*cos(lon);
y = cos(lat).*sin(lon);
z = sin(lat);

v = [x(:), y(:), z(:)];
end