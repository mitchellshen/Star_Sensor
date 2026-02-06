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
d_fld = 'imap_lo/';
% CDF  folder
c_fld = 'matlab_cdf391_patch-arm64/';

addpath([dir d_fld]);
addpath([dir c_fld]);

%==============================================
%% Import & Process CDFs
close all;
CDFdate = ['20251226';'20251227';'20251228';'20251229';'20251230';...
           '20260121';'20260129';'20260130'];

CDFrept = [  90, 91, 92, 93, 94,...
            133,141,142];

CDFpivt = [  90, 75, 90,105, 90,...
             90, 90, 90];

for i = 7%1:size(CDFdate,1)

idx = i;

[l1a_star_data, l1a_star_info] = ...
    spdfcdfread(['imap_lo_l1a_star_' CDFdate(idx,:) ...
                 '-repoint' num2str(CDFrept(idx),'%05d') '_v002.cdf']);

l1a_star_var = l1a_star_info.Variables(:, [1 4]);


[l1b_star_data, l1b_star_info] = ...
    spdfcdfread(['imap_lo_l1b_prostar_' CDFdate(idx,:) ...
                 '-repoint' num2str(CDFrept(idx),'%05d') '_v002.cdf']);

l1b_star_var = l1b_star_info.Variables(:, [1 4]);

l1a_star_shcoarse   = cell2mat(l1a_star_data( 8));
l1a_star_count      = cell2mat(l1a_star_data( 9));
l1a_star_chksum     = cell2mat(l1a_star_data(10));
l1a_star_datadata   = cell2mat(l1a_star_data(11));
l1a_star_datadata(abs(l1a_star_datadata) > 1e30) = NaN; % fill values
l1a_star_epoch      = cell2mat(l1a_star_data(12));
l1a_star_data_index = cell2mat(l1a_star_data(13));

l1b_star_spin_ang_bin = cell2mat(l1b_star_data(1));
l1b_star_met          = cell2mat(l1b_star_data(2));
l1b_star_data_ave_amp = cell2mat(l1b_star_data(3));
l1b_star_data_ave_amp(abs(l1b_star_data_ave_amp) > 1e30) = NaN; % fill values
l1b_star_cunt_per_bin = cell2mat(l1b_star_data(4));
l1b_star_epoch        = cell2mat(l1b_star_data(5));
l1b_star_data_spn_ang = cell2mat(l1b_star_data(6));


%==============================================
%% Ave. STAR oplot
figure(1)

plot(l1b_star_data_spn_ang, l1b_star_data_ave_amp); hold on;

plot(l1b_star_data_spn_ang, mean(l1b_star_data_ave_amp), 'k-', 'LineWidth', 1.6); 
hold on;

xlabel('NEP Angle (degree)');
ylabel('Star Sensor Amplitude (mV)');
xlim([0 360]);
ylim([0 2500]);
labels = [arrayfun(@(k) sprintf('spin_ave_#%d',k),...
    1:size(l1b_star_data_ave_amp,1), 'UniformOutput', false),...
    {'Ave of L1b star'}];
lgd = legend(labels,'Location','eastoutside', 'Interpreter','none');
lgd.NumColumns = 2;
set(gca,'lineWidth',1,'FontSize',12,'FontWeight','default');
grid on; box on;

%% plot individual day
figure(2)

% L1a
% plot(double(mod(l1a_star_data_index+124,720))/2, mean(l1a_star_datadata), ...
%      'LineStyle','none', 'Marker','+', 'MarkerSize',5);
plot(mod(l1b_star_data_spn_ang+62,360), mean(l1a_star_datadata), ...
     'LineStyle','none', 'Marker','+', 'MarkerSize',5); hold on;

% L1b
plot(l1b_star_data_spn_ang, mean(l1b_star_data_ave_amp), ...
     'LineStyle','none', 'Marker','x', 'MarkerSize',5); hold on;

grid on;
xlabel('NEP Angle (degree)');
ylabel('Star Sensor Amplitude (mV)');
xlim([0 360]);
ylim([0 2500]);
% label2{idx} = [CDFdate(idx,:) ', pivot ' num2str(CDFpivt(idx),'%03d')];
% lgd = legend(label2,'Location','best', 'Interpreter','none');
% lgd.NumColumns = 2;
legend({'l1b','l1a'},'Location','best', 'Interpreter','none');
set(gca,'lineWidth',1,'FontSize',12,'FontWeight','default');
grid on; box on;

% plot multiple days comparison
figure(3)

% L1b
plot(l1b_star_data_spn_ang, mean(l1b_star_data_ave_amp), ...
     'LineStyle','none', 'Marker','x', 'MarkerSize',5); hold on;

grid on;
xlabel('NEP Angle (degree)');
ylabel('Star Sensor Amplitude (mV)');
xlim([0 360]);
ylim([0 2500]);
label2{idx} = [CDFdate(idx,:) ', pivot ' num2str(CDFpivt(idx),'%03d')];
lgd = legend(label2,'Location','best', 'Interpreter','none');
lgd.NumColumns = 2;
set(gca,'lineWidth',1,'FontSize',12,'FontWeight','default');
grid on; box on;

% end


%% 2D colormap 

t = datetime(l1b_star_epoch,'ConvertFrom','datenum','TimeZone','UTC');
Z = l1b_star_data_ave_amp.';   % (spin angle × time)

% continued
figure(5)
h = imagesc(t, l1b_star_data_spn_ang, Z);
axis xy
ylim([0 360]);
hcb = colorbar;
clim([0 2500]);
colormap(parula);
xlabel('Epoch');
ylabel(['NEP Angle (degree)']);
ylabel(hcb,'Star Sensor Amplitude (mV)','FontSize',12);

% Make NaNs transparent (NOT colored)
set(h,'AlphaData', ~isnan(Z))

% Set axes background to white → NaNs appear white
set(gca,'Color','w')

set(gca,'lineWidth',1,'FontSize',12,'FontWeight','default');


% descrete
figure(6)
[T, ANG] = ndgrid(l1b_star_epoch, l1b_star_data_spn_ang);
x = T(:);                          % epoch
y = ANG(:);                        % spin angle
c = l1b_star_data_ave_amp(:);      % amplitude
scatter(datetime(x,'ConvertFrom','datenum'), y, 20, c, 'filled')
xlabel('Epoch')
ylabel('NEP Angle (degree)')
cb = colorbar;
ylabel(cb,'Star Sensor Amplitude (mV)','FontSize',12)
box on; grid on;
ylim([0 360]);
set(gca,'lineWidth',1,'FontSize',12,'FontWeight','default');

end


%% Autocorrelation test
% not necessary good checks due to detrending
% 20251226, pivot 090Best lag = -58.08 deg (corr=0.958)
% 20251227, pivot 075Best lag = -59.07 deg (corr=0.958)
% 20251228, pivot 090Best lag = -58.07 deg (corr=0.958)
% 20251229, pivot 105Best lag = -62.57 deg (corr=0.971)
% 20251230, pivot 090Best lag = -58.07 deg (corr=0.958)
% 20260121, pivot 090Best lag = -63.17 deg (corr=0.937)

% Angle grid (degrees). If your l1b_star_data_spn_ang already is 0:0.5:359.5, use it.
ang = l1b_star_data_spn_ang(:);

y_l1b = mean(l1b_star_data_ave_amp, 1, 'omitnan');   % 1×720 (if ave_amp is 31×720)
y_l1b = y_l1b(:);

y_l1a = mean(l1a_star_datadata, 1, 'omitnan');       % should be 1×720
y_l1a = y_l1a(:);

% If either has NaNs, fill them to avoid xcorr choking.
% (Simple linear fill; ok for correlation diagnostics)
y_l1a = fillmissing(y_l1a,'linear','EndValues','nearest');
y_l1b = fillmissing(y_l1b,'linear','EndValues','nearest');

% Remove mean + normalize (so correlation compares SHAPE, not amplitude scale)
y_l1a = (y_l1a - mean(y_l1a)) / std(y_l1a);
y_l1b = (y_l1b - mean(y_l1b)) / std(y_l1b);

dx = median(diff(ang));   % degrees per sample, likely 0.5


[ac_a, lags] = xcorr(y_l1a, 'coeff');
[ac_b, ~]    = xcorr(y_l1b, 'coeff');
lag_deg = lags * dx;

% Plot (focus on +/-180 deg)
use = (lag_deg >= -180) & (lag_deg <= 180);

figure(8)
plot(lag_deg(use), ac_a(use), 'LineWidth', 1.5); hold on
plot(lag_deg(use), ac_b(use), 'LineWidth', 1.5); hold on
grid on
xlabel('Lag (deg)')
ylabel('Autocorrelation (coeff)')
title('Autocorrelation: L1A vs L1B')
legend({'L1A','L1B'}, 'Location','eastoutside', 'Interpreter','none')


[cc, lags] = xcorr(y_l1a, y_l1b, 'coeff');  % L1A vs L1B
lag_deg = lags * dx;

use = (lag_deg >= -180) & (lag_deg <= 180);

figure(9)
plot(lag_deg(use), cc(use), 'LineWidth', 1.5); hold on
grid on
xlabel('Lag (deg)')
ylabel('Cross-correlation (coeff)')
title('Cross-correlation: L1A vs L1B')


[ccmax, imax] = max(cc);
best_lag_deg  = lag_deg(imax);

fprintf([CDFdate(idx,:) ', pivot ' num2str(CDFpivt(idx),'%03d')...
    'Best lag = %.2f deg (corr=%.3f)\n'], best_lag_deg, ccmax);

