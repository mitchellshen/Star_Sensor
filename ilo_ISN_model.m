%% MAIN
close all; clearvars; clc; format long; warning('off','all')

%% Root folder
fld = ['/Users/mitchellshen/Library/CloudStorage/Dropbox-Personal/' ...
       'IMAP_ILO/SCI_model/' ...
       '20260421_Modeling of ISN update of ionization/' ...
       'ion.2026c.patched/'];

assert(exist(fld, 'dir') == 7, 'Root folder does not exist.')

outdir = [fld 'plots'];
if ~exist(outdir,'dir')
    mkdir(outdir);
end

%% Find all .dat files recursively
files = dir(fullfile(fld, '**', '*.dat'));
fprintf('Found %d .dat files\n', numel(files));

%% Initialize nested output structure
sp_list  = {'He','H','O','Ne','unknown'};
pop_list = {'pr','sc','unknown'};

all_data = struct();
for is = 1:numel(sp_list)
    for ip = 1:numel(pop_list)
        all_data.(sp_list{is}).(pop_list{ip}) = {};
    end
end

%% Preallocate metadata containers
N = numel(files);

species_col = strings(N,1);
pop_col     = strings(N,1);
name_col    = strings(N,1);
folder_col  = strings(N,1);
path_col    = strings(N,1);

repoint_col = nan(N,1);
timeyr_col  = nan(N,1);
angle_col   = nan(N,1);

nrow_col    = nan(N,1);
ncol_col    = nan(N,1);

read_ok_col = false(N,1);
msg_col     = strings(N,1);

%% Read all files and build metadata
for i = 1:N
    file_fullpath = fullfile(files(i).folder, files(i).name);
    folder_lower  = lower(files(i).folder);
    name_lower    = lower(files(i).name);
    parts         = split(folder_lower, filesep);

    %-------------------------
    % species classification
    %-------------------------
    if contains(folder_lower, 'helium')
        species = "He";
    elseif contains(folder_lower, 'hydrogen')
        species = "H";
    elseif contains(folder_lower, 'oxygen')
        species = "O";
    elseif contains(folder_lower, 'neon')
        species = "Ne";
    else
        species = "unknown";
    end

    %-------------------------
    % population classification
    %-------------------------
    if any(strcmp(parts, 'pr')) || contains(name_lower, '_pri_')
        pop = "pr";
    elseif any(strcmp(parts, 'sc')) || contains(name_lower, '_sec_')
        pop = "sc";
    else
        pop = "unknown";
    end

    %-------------------------
    % store basic metadata
    %-------------------------
    species_col(i) = species;
    pop_col(i)     = pop;
    name_col(i)    = string(files(i).name);
    folder_col(i)  = string(files(i).folder);
    path_col(i)    = string(file_fullpath);

    %-------------------------
    % parse repoint / time / angle from filename
    % Example:
    % simColl6d_He_pri_grid_repoint00077_2025.950_105.0.dat
    %-------------------------
    tok = regexp(files(i).name, ...
        'repoint0*(\d+)_([0-9]+\.[0-9]+)_([0-9]+\.[0-9]+)\.dat$', ...
        'tokens', 'once');

    if ~isempty(tok)
        repoint_col(i) = str2double(tok{1});
        timeyr_col(i)  = str2double(tok{2});
        angle_col(i)   = str2double(tok{3});
    end

    fprintf('[%4d/%4d] %s | %s | %s\n', i, N, species, pop, files(i).name);

    %-------------------------
    % read actual spin-table
    %-------------------------
    try
        data = read_isn_dat(file_fullpath);

        entry = struct();
        entry.name    = files(i).name;
        entry.folder  = files(i).folder;
        entry.path    = file_fullpath;
        entry.species = char(species);
        entry.pop     = char(pop);
        entry.repoint = repoint_col(i);
        entry.time_yr = timeyr_col(i);
        entry.angle   = angle_col(i);
        entry.data    = data;

        all_data.(char(species)).(char(pop)){end+1} = entry;

        nrow_col(i) = size(data,1);
        ncol_col(i) = size(data,2);

        read_ok_col(i) = true;
        msg_col(i)     = "OK";

    catch ME
        read_ok_col(i) = false;
        msg_col(i)     = string(ME.message);

        fprintf('    FAILED: %s\n', ME.message);
    end
end

%% Build metadata table
T = table(species_col, pop_col, name_col, folder_col, path_col, ...
          repoint_col, timeyr_col, angle_col, ...
          nrow_col, ncol_col, read_ok_col, msg_col, ...
    'VariableNames', {'species','pop','name','folder','path', ...
                      'repoint','time_yr','angle', ...
                      'nrow','ncol','read_ok','message'});

T = sortrows(T, {'species','pop','repoint','time_yr','angle'});

%% Summary
fprintf('\n===== SUMMARY =====\n');
fprintf('Readable files: %d / %d\n', sum(T.read_ok), height(T));

disp(groupsummary(T, {'species','pop'}, 'sum', 'read_ok'))

disp(unique(T(:,{'nrow','ncol'}), 'rows'))

%% Optional save
% writetable(T, 'isn_model_file_index.csv')
% save('isn_model_all_data.mat', 'all_data', 'T', '-v7.3')

%% Modeling Plot
% ==========================
% Plot one repoint in 2x3 panels
% Top row: flux vs spin angle
% Bottom row: energy [eV] vs spin angle
% ==========================
close all;
rep_list = sort(unique(T.repoint(T.read_ok))); % Loop over all repoints
rep_list = rep_list(~isnan(rep_list));

pivot_targets = [75 90 105];   % always show these 3 panels
tol = 3;                       % allow nearby pivots, e.g. 76 for 75, 104 for 105

yr = floor(T.time_yr);
frac = T.time_yr - yr;
isLeap = (mod(yr,4)==0 & mod(yr,100)~=0) | (mod(yr,400)==0);
daysInYear = 365 + isLeap;
doy = floor(frac .* daysInYear) + 1;
doy(doy > daysInYear) = daysInYear(doy > daysInYear);

% column definitions from read_isn_dat output:
% 1 = Spin
% 2 = flux
% 3 = speed
% 4 = speed^2
% 5 = speed^3
% 6 = long_ecl
% 7 = lat_ecl
spin_col   = 1;
flux_col   = 2;
speed2_col = 4;


% physical constants
amu = 1.66053906660e-27;    % kg
qe  = 1.602176634e-19;      % J/eV

% species masses in atomic mass units
mass_amu.He = 4;
mass_amu.H  = 1;
mass_amu.O  = 16;
mass_amu.Ne = 20;

% color setup
C.He_pr = [0.85 0.45 0.10];   % orange
C.He_sc = [1.00 0.00 1.00];   % magenta
C.H_pr  = [0.10 0.10 0.10];   % black
C.H_sc  = [0.60 0.60 0.60];   % gray
C.O_pr  = [0.00 0.55 0.45];   % teal
C.O_sc  = [0.20 0.90 0.20];   % green
C.Ne_pr = [0.20 0.70 1.00];   % cyan
C.unk   = [0.50 0.50 0.50];

for ir = 1:numel(rep_list)

    rep0 = rep_list(ir);
    fprintf('\n==============================\n');
    fprintf('Repoint %d (%d/%d)\n', rep0, ir, numel(rep_list));
    fprintf('==============================\n');

    % ----- repoint-level date string -----
    idx_rep = T.read_ok & T.repoint == rep0;
    
    yr_rep  = unique(yr(idx_rep));
    doy_rep = unique(doy(idx_rep));
    
    yr_rep  = yr_rep(~isnan(yr_rep));
    doy_rep = doy_rep(~isnan(doy_rep));
    
    if ~isempty(yr_rep) && ~isempty(doy_rep)
        date_str = sprintf('%d-%03d', yr_rep(1), doy_rep(1));
    else
        date_str = 'NaN';
    end


    figure('Visible','off','Color','w','Position',[100 100 1100 650]);
    tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

    for ia = 1:numel(pivot_targets)
        piv0 = pivot_targets(ia);

        % top axis = flux
        ax1 = nexttile(ia);
        hold(ax1,'on')
        grid(ax1,'on')
        box(ax1,'on')
        set(ax1,'YScale','log')

        % bottom axis = energy
        ax2 = nexttile(ia+3);
        hold(ax2,'on')
        grid(ax2,'on')
        box(ax2,'on')
        set(ax2,'YScale','log')

        % match nearby pivot values
        idx = T.read_ok & T.repoint == rep0 & abs(T.angle - piv0) <= tol;
        Ta = T(idx,:);
        

        fprintf('\n--- repoint %d, target pivot %.1f ---\n', rep0, piv0);
        fprintf('Matched files: %d\n', height(Ta));

        labels_plotted = {};

        for i = 1:height(Ta)
            sp  = char(Ta.species(i));
            pop = char(Ta.pop(i));
            pth = char(Ta.path(i));

            data = read_isn_dat(pth);

            if size(data,2) < max([spin_col flux_col speed2_col])
                continue
            end

            x      = data(:,spin_col);
            y_flux = data(:,flux_col);
            v2     = data(:,speed2_col);

            % convert speed^2 -> energy [eV]
            % assumes speed is in km/s, so convert (km/s)^2 -> (m/s)^2
            if isfield(mass_amu, sp)
                m = mass_amu.(sp) * amu;
            else
                m = NaN;
            end

            y_E = 0.5 * m .* (v2 * 1e6) / qe;   % eV

            good_flux = isfinite(x) & isfinite(y_flux) & (y_flux > 0);
            good_E    = isfinite(x) & isfinite(y_E)    & (y_E > 0);

            if ~any(good_flux) && ~any(good_E)
                continue
            end

            % label + color
            if strcmp(sp,'He') && strcmp(pop,'pr')
                clr = C.He_pr; lbl = 'He_{pri}';
            elseif strcmp(sp,'He') && strcmp(pop,'sc')
                clr = C.He_sc; lbl = 'He_{sec}';
            elseif strcmp(sp,'H') && strcmp(pop,'pr')
                clr = C.H_pr;  lbl = 'H_{pri}';
            elseif strcmp(sp,'H') && strcmp(pop,'sc')
                clr = C.H_sc;  lbl = 'H_{sec}';
            elseif strcmp(sp,'O') && strcmp(pop,'pr')
                clr = C.O_pr;  lbl = 'O_{pri}';
            elseif strcmp(sp,'O') && strcmp(pop,'sc')
                clr = C.O_sc;  lbl = 'O_{sec}';
            elseif strcmp(sp,'Ne')
                clr = C.Ne_pr; lbl = 'Ne_{pri}';
            else
                clr = C.unk;   lbl = [sp '_' pop];
            end

            % top plot
            if any(good_flux)
                if ~ismember(lbl, labels_plotted)
                    plot(ax1, x(good_flux), y_flux(good_flux), ...
                        'LineWidth', 1.8, 'Color', clr, 'DisplayName', lbl);
                    labels_plotted{end+1} = lbl;
                else
                    plot(ax1, x(good_flux), y_flux(good_flux), ...
                        'LineWidth', 1.8, 'Color', clr, 'HandleVisibility','off');
                end
            end

            % bottom plot
            if any(good_E)
                plot(ax2, x(good_E), y_E(good_E), ...
                    'LineWidth', 1.8, 'Color', clr, 'HandleVisibility','off');
            end
        end

        % formatting
        xlim(ax1,[0 360])
        xlim(ax2,[0 360])
        ylim(ax1,[1e-1 1e7])
        ylim(ax2,[1e0  1e3])
        ylabel(ax1,'flux')
        ylabel(ax2,'energy [eV]')
        xlabel(ax2,'spin angle')

        title(ax1, sprintf('DOY = %s, repoint = %d, \\epsilon = %.0f', ...
            date_str, rep0, piv0))

        if ~isempty(labels_plotted)
            legend(ax1,'Location','best')
        else
            text(ax1, 0.5, 0.5, 'no data', ...
                'Units','normalized', ...
                'HorizontalAlignment','center', ...
                'FontSize',12)
            text(ax2, 0.5, 0.5, 'no data', ...
                'Units','normalized', ...
                'HorizontalAlignment','center', ...
                'FontSize',12)
        end
    end

    drawnow

    % ===== SAVE FIGURE =====
    tag = 'ISN_model';
    
    fname_png = fullfile(outdir, ...
        sprintf('%s_%s_repoint_%05d.png', tag, date_str, rep0));
    
    fname_pdf = fullfile(outdir, ...
        sprintf('%s_%s_repoint_%05d.pdf', tag, date_str, rep0));
    
    % exportgraphics(gcf, fname_png, 'Resolution',300);
    % exportgraphics(gcf, fname_pdf);
    % ===== CLOSE FIGURE (recommended) =====
    close(gcf)

end


%% Sky map
%%% Final accumulated pseudo skymaps by pivot angle
% Sum all species/populations together
% Combine over ALL repoints, but KEEP pivot groups separate

close all;

outdir = [fld 'pseudo_skymap_pivot'];
if ~exist(outdir,'dir')
    mkdir(outdir);
end

% ----- pivot groups -----
pivot_targets = [75 90 105];
tol = 3;   % e.g. 76 goes into 75-bin, 104 into 105-bin

% ----- sky grid -----
lon_edges = 0:2:360;
lat_edges = -90:2:90;

lon_cent = lon_edges(1:end-1) + diff(lon_edges)/2;
lat_cent = lat_edges(1:end-1) + diff(lat_edges)/2;

% store final maps
map_final = nan(numel(lat_cent), numel(lon_cent), numel(pivot_targets));

for ip = 1:numel(pivot_targets)

    piv0 = pivot_targets(ip);

    % all files across all repoints within this pivot group
    idx = T.read_ok & abs(T.angle - piv0) <= tol;
    Ta = T(idx,:);

    fprintf('\n==============================\n');
    fprintf('Accumulating all repoints for pivot ~%.0f\n', piv0);
    fprintf('Matched files: %d\n', height(Ta));
    fprintf('==============================\n');

    all_lon  = [];
    all_lat  = [];
    all_flux = [];

    for i = 1:height(Ta)
        pth = char(Ta.path(i));
        data = read_isn_dat(pth);

        if size(data,2) < 7
            continue
        end

        flux = data(:,2);
        lon  = data(:,6);
        lat  = data(:,7);

        good = isfinite(flux) & isfinite(lon) & isfinite(lat) & (flux > 0);

        all_flux = [all_flux; flux(good)];
        all_lon  = [all_lon;  lon(good)];
        all_lat  = [all_lat;  lat(good)];
    end

    if isempty(all_flux)
        fprintf('No valid data for pivot ~%.0f\n', piv0);
        continue
    end

    % bin onto sky grid
    [~,~,ix] = histcounts(all_lon, lon_edges);
    [~,~,iy] = histcounts(all_lat, lat_edges);

    goodbin = ix > 0 & iy > 0;

    map_sum = accumarray([iy(goodbin), ix(goodbin)], all_flux(goodbin), ...
        [numel(lat_cent), numel(lon_cent)], @sum, NaN);

    map_final(:,:,ip) = map_sum;
end

% ----------------------------
% Figure 1: three final maps
% ----------------------------
figure('Color','w','Position',[50 50 1400 300]);
% 'Visible','off',
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

lon_center = -105;

for ip = 1:numel(pivot_targets)
    ax = nexttile;
    piv0 = pivot_targets(ip);
    thismap = map_final(:,:,ip);

    % --- recenter longitude around lon_center ---
    lon_plot = mod(lon_cent - lon_center + 180, 360) - 180 + lon_center;

    % --- sort longitude and reorder map accordingly ---
    [lon_plot, idx_sort] = sort(lon_plot);
    thismap_plot = thismap(:, idx_sort);

    imagesc(lon_plot, lat_cent, thismap_plot);
    set(ax,'YDir','normal','XDir','reverse');
    set(ax,'ColorScale','log');
    colormap(ax, turbo);


    xlabel(ax,'long\_ecl [deg]');
    ylabel(ax,'lat\_ecl [deg]');
    title(ax, sprintf('Combined over all repoints, pivot \\approx %.0f', piv0));
    grid(ax,'on');
    ax.GridColor = [0.7 0.7 0.7]; ax.GridAlpha = 0.25;

    % optional but useful so the center is really at lon_center
    xlim(ax,[lon_center-180, lon_center+180])
    ylim(ax,[-90  90])
    clim(ax,[1e1 1e6])
    
    % --- custom x ticks and wrapped labels ---
    xt = lon_center + (-180:60:180);
    xticks(ax, xt)
    xtlbl = mod(xt + 180, 360) - 180;
    xticklabels(ax, string(round(xtlbl)))
    
    % --- display y ticks in usual way ---
    yticks(ax, -90:30:90);

    cb = colorbar(ax);
    cb.Label.String = 'summed flux';
    

    set(ax,'Color',[0.2 0.2 0.2]);
    h = findobj(ax,'Type','Image');
    if ~isempty(h)
        set(h, 'AlphaData', ~isnan(thismap_plot))
    end

    if all(isnan(thismap_plot(:)))
        text(ax, 0.5, 0.5, 'no data', ...
            'Units','normalized', ...
            'HorizontalAlignment','center', ...
            'FontSize',14)
    end
end

fname_png = fullfile(outdir, 'final_accumulated_maps_by_pivot.png');
fname_pdf = fullfile(outdir, 'final_accumulated_maps_by_pivot.pdf');

% exportgraphics(gcf, fname_png, 'Resolution', 300);
% exportgraphics(gcf, fname_pdf);
% close(gcf)

% ----------------------------
% Optional: save each pivot map separately too
% ----------------------------
for ip = 1:numel(pivot_targets)
    piv0 = pivot_targets(ip);
    thismap = map_final(:,:,ip);

    figure('Visible','off','Color','w','Position',[100 100 780 480]);

    imagesc(lon_cent, lat_cent, thismap);
    set(gca,'YDir','normal');
    set(gca,'ColorScale','log');
    colormap(turbo);

    xlabel('long\_ecl [deg]');
    ylabel('lat\_ecl [deg]');
    title(sprintf('Combined over all repoints, pivot \\approx %.0f', piv0));
    grid on

    cb = colorbar;
    cb.Label.String = 'summed flux';

    set(gca,'Color',[1 1 1]);
    h = findobj(gca,'Type','Image');
    if ~isempty(h)
        set(h, 'AlphaData', ~isnan(thismap))
    end

    if all(isnan(thismap(:)))
        text(0.5, 0.5, 'no data', ...
            'Units','normalized', ...
            'HorizontalAlignment','center', ...
            'FontSize',14)
    end

    fname_png = fullfile(outdir, sprintf('final_pivot_%03d.png', round(piv0)));
    fname_pdf = fullfile(outdir, sprintf('final_pivot_%03d.pdf', round(piv0)));
    % 
    % exportgraphics(gcf, fname_png, 'Resolution', 300);
    % exportgraphics(gcf, fname_pdf);
    % close(gcf)
end

%% Mollweide map (NEED REVISE)

% ----------------------------
% Figure 2: three final maps in Mollweide projection
% ----------------------------
figure('Color','w','Position',[50 380 1400 420]);
% 'Visible','off',

lon_center = -105;
lat_center = 5;

for ip = 1:numel(pivot_targets)

    subplot(1,3,ip)

    piv0 = pivot_targets(ip);
    thismap = map_final(:,:,ip);

    % --- recenter longitude around lon_center ---
    lon_plot = mod(lon_cent - lon_center + 180, 360) - 180 + lon_center;

    % --- sort longitude and reorder map accordingly ---
    [lon_plot, idx_sort] = sort(lon_plot);
    thismap_plot = thismap(:, idx_sort);

    % --- build lon/lat mesh for Mapping Toolbox ---
    [LonGrid, LatGrid] = meshgrid(lon_plot, lat_cent);

    % --- create Mollweide map axes ---
    axesm('mollweid', ...
    'Origin', [lat_center lon_center 0], ...
    'MapLatLimit', [-90 90], ...
    'Frame', 'on', ...
    'Grid', 'on', ...
    'MeridianLabel', 'off', ...
    'ParallelLabel', 'off', ...
    'MLineLocation', 30, ...
    'PLineLocation', 30, ...
    'GLineStyle', ':', ...
    'Gcolor', [0.5 0.5 0.5], ...
    'GLineWidth', 0.6, ...
    'FEdgeColor', [0.2 0.2 0.2], ...
    'FLineWidth', 1.0);

    % --- plot gridded map ---
    h = surfm(LatGrid, LonGrid, thismap_plot);
    tightmap
    set(gca, 'FontSize', 11)

    % make NaNs transparent / background white
    set(gca, 'Color', [1 1 1]);
    if ~isempty(h)
        set(h, 'AlphaData', ~isnan(thismap_plot));
    end

    % colormap / color scale
    colormap(gca, turbo);
    set(gca, 'ColorScale', 'log');

    title(sprintf('Combined over all repoints, pivot \\approx %.0f', piv0), ...
        'FontWeight', 'bold');

    cb = colorbar;
    cb.Label.String = 'summed flux';

    if all(isnan(thismap_plot(:)))
        text(0.5, 0.5, 'no data', ...
            'Units','normalized', ...
            'HorizontalAlignment','center', ...
            'FontSize',14)
    end
end

fname_png = fullfile(outdir, 'final_accumulated_maps_by_pivot_mollweide.png');
fname_pdf = fullfile(outdir, 'final_accumulated_maps_by_pivot_mollweide.pdf');

% exportgraphics(gcf, fname_png, 'Resolution', 300);
% exportgraphics(gcf, fname_pdf);
% close(gcf)



%% FUNCTIOn
function data = read_isn_dat(file_fullpath)

    fid = fopen(file_fullpath, 'r');
    if fid < 0
        error('Cannot open file: %s', file_fullpath);
    end

    cleanup = onCleanup(@() fclose(fid));

    % Find the line that contains the real table header
    found_header = false;
    while ~feof(fid)
        thisline = fgetl(fid);
        if ischar(thisline) && contains(thisline, '#Spin')
            found_header = true;
            break
        end
    end

    if ~found_header
        error('Did not find #Spin header in file: %s', file_fullpath);
    end

    % Read the numeric table after the #Spin line
    C = textscan(fid, '%f %f %f %f %f %f %f', ...
        'MultipleDelimsAsOne', true, ...
        'CollectOutput', true);

    data = C{1};

    if isempty(data)
        error('No numeric table found after #Spin in file: %s', file_fullpath);
    end
end