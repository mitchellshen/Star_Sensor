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
%% Import CDF & CSV
[l1b_star_data, l1b_star_info] = ...
    spdfcdfread('imap_lo_l1b_star_20251110_v999.cdf');
    % spdfcdfread('imap_lo_l1b_instrument-status-summary_20251104-repoint00038_v001.cdf');

l1b_star_data_spn_ang = cell2mat(l1b_star_data(6));
l1b_star_data_ave_amp = cell2mat(l1b_star_data(3));
l1b_star_data_ave_amp(abs(l1b_star_data_ave_amp) > 1e30) = NaN;
l1b_star_spin_ang_bin = cell2mat(l1b_star_data(1));
l1b_star_cunt_per_bin = cell2mat(l1b_star_data(4));
l1b_star_met          = cell2mat(l1b_star_data(2));
l1b_star_epoch        = cell2mat(l1b_star_data(5));

avg_star = readtable('imap_lo_avg_star_2025314.csv');
avg_star_deg = table2array(avg_star(3:end,1));
avg_star_amp = table2array(avg_star(3:end,2));

% l1b_star_var = l1b_star_info.Variables;
l1b_star_var = l1b_star_info.Variables(:, [1 4]);




%==============================================
%% Ave. STAR oplot
figure(1)

plot(l1b_star_data_spn_ang, l1b_star_data_ave_amp); hold on;

plot(avg_star_deg, avg_star_amp, 'k-.', 'LineWidth', 2.0)
xlabel('Degree');
ylabel('Amplitude');
title('Average STAR Data');
xlim([0 360])
grid on; hold on;

plot(l1b_star_data_spn_ang, mean(l1b_star_data_ave_amp), 'b-.', 'LineWidth', 2.0); 
hold on;

labels = [arrayfun(@(k) sprintf('spin_ave_#%d',k),...
    1:size(l1b_star_data_ave_amp,1), 'UniformOutput', false),...
    {'Raw CSV'}, {'Ave of L1b star'}];
legend(labels,'Location','eastoutside', 'Interpreter','none');


figure(2)

plot(avg_star_deg, avg_star_amp, 'k-.', 'LineWidth', 2.0)
xlabel('Degree');
ylabel('Amplitude');
title('Average STAR Data');
xlim([0 360])
grid on; hold on;

plot(l1b_star_data_spn_ang, mean(l1b_star_data_ave_amp), 'b-.', 'LineWidth', 2.0); 
hold on;

labels = [{'Raw CSV'}, {'Ave of L1b star'}];
legend(labels,'Location','eastoutside', 'Interpreter','none');

%% Ave. STAR auto-corr
%  Raw vs. Ave of (64 spin-ave.) = 31 rows data

% --- Inputs (edit these to your variable names) ---
x_raw = avg_star_deg;      % degrees for raw CSV curve
y_raw = avg_star_amp;      % raw CSV amplitude (black)

x_l1b = l1b_star_data_spn_ang;               % degrees for L1b averaged curve
y_l1b = l1b_star_data_ave_amp;               % L1b averaged amplitude (blue)

% --- Make a common uniform grid in degrees ---
dx = 1;                              % choose 1 deg (or your native step)
xg = (0:dx:360-dx).';                % column vector, 0..359

% Ensure degrees wrap consistently (optional but usually helpful)
x_raw = mod(x_raw,360);
x_l1b = mod(x_l1b,360);

% Sort (interp1 requires monotonic x)
[x_raw, i1] = sort(x_raw);  y_raw = y_raw(i1);
[x_l1b, i2] = sort(x_l1b);  y_l1b = y_l1b(i2);

% Interpolate onto common grid
yr = interp1(x_raw, y_raw, xg, 'linear', 'extrap');
yb = interp1(x_l1b, y_l1b, xg, 'linear', 'extrap');

% --- Clean / normalize (so comparison is fair) ---
yr = yr - mean(yr,'omitnan');
yb = yb - mean(yb,'omitnan');

yr = yr ./ std(yr,'omitnan');
yb = yb ./ std(yb,'omitnan');

% --- Autocorrelation ---
[r_raw, lags] = xcorr(yr, 'coeff');
[r_l1b, ~   ] = xcorr(yb, 'coeff');
lag_deg = lags * dx;

% Optional: only show +-180 deg to focus on structure
use = (lag_deg >= -180) & (lag_deg <= 180);

% --- Plot ---
figure(3)
plot(lag_deg(use), r_raw(use), 'k--', 'LineWidth', 2); hold on
plot(lag_deg(use), r_l1b(use), 'b-.', 'LineWidth', 2);
grid on
xlabel('Lag (degree)')
ylabel('Autocorrelation (coeff)')
title('Autocorrelation: Raw CSV vs Spin-averaged L1b')
legend({'Raw CSV','Spin-averaged L1b'}, 'Location','eastoutside', 'Interpreter','none')

% --- Optional: difference curve (how much smoothing changed structure) ---
figure(4)
plot(lag_deg(use), (r_l1b(use) - r_raw(use)), 'k', 'LineWidth', 2);
grid on
xlabel('Lag (degree)')
ylabel('\Delta autocorr (L1b - Raw)')
title('Autocorrelation difference')

%% 2D colormap
t = datetime(l1b_star_epoch,'ConvertFrom','datenum','TimeZone','UTC');
Z = l1b_star_data_ave_amp.';   % (spin angle × time)

figure(5)
h = imagesc(t, l1b_star_data_spn_ang, Z);
axis xy
colorbar
colormap(parula)     % or 'parula', etc.
xlabel('Epoch')

% Make NaNs transparent (NOT colored)
set(h,'AlphaData', ~isnan(Z))

% Set axes background to white → NaNs appear white
set(gca,'Color','w')


figure(6)
[T, ANG] = ndgrid(l1b_star_epoch, l1b_star_data_spn_ang);
x = T(:);                          % epoch
y = ANG(:);                        % spin angle
c = l1b_star_data_ave_amp(:);      % amplitude
scatter(datetime(x,'ConvertFrom','datenum'), y, 20, c, 'filled')
xlabel('Epoch time')
ylabel('Spin angle (deg)')
colorbar


