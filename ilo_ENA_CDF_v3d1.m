%% MAIN

close all; clearvars; clc; format long; warning('off','all')

fld = '/Users/mitchellshen/CAVA_IDL/IMAP-Lo June 10th 2026 Maps/l3/2025/11/';
cdf_name = 'imap_lo_l3_ilo-enasbsMsk-h-hf-sp-ram-hae-6deg-6mo_20251125_v001.cdf';
cdf_file = [fld cdf_name];

%==================== USER SWITCH ====================%
flip_view = true;   % true = outside-looking-in style
slice_view= false;
rate      = false;  % kept only for downstream colorbar/clim branching (see below)
save      = false;
esa_step  = 4;       % index into the file's 'energy' variable (1..7 for this product)
epoch_idx = 1;       % which time record to plot, if the CDF has more than one
quantity  = 'flux';  % 'flux' -> ena_intensity, 'bg' -> bg_intensity, 'survival' -> survival_probability
esa_volt  = {'16','30','56','106','200','404','787'};  % fallback labels only; real energy is read from the file below

%==================== READ DATA (FROM CDF) ====================%
% This CDF (a 6-month IMAP-Lo survival-corrected ENA map) stores:
%   longitude, latitude   -> pixel-CENTER coordinates (deg), irregular counts
%                             (e.g. this file has 58 longitude x 30 latitude
%                             bins -- NOT a fixed 60x30 grid like the CSV maps)
%   energy                -> 7 energy passband centers (keV), DEPEND_1 of ena_intensity
%   ena_intensity, bg_intensity, survival_probability, ena_intensity_stat_uncert, ...
%                          -> DEPEND_0..3 = [epoch, energy, longitude, latitude]
%
% cdfread's returned array dimension order for a given MATLAB/CDF library
% version isn't something I can verify without running MATLAB myself, so
% instead of hard-coding [energy, lon, lat] we MATCH each returned
% dimension's length against numel(energy)/numel(longitude)/numel(latitude).
% Those three are 7 / 58 / 30 here -- all different -- so the match is
% unambiguous and self-correcting even if the library's convention changes.

switch quantity
    case 'flux',     zvar = 'ena_intensity';
    case 'bg',        zvar = 'bg_intensity';
    case 'survival', zvar = 'survival_probability';
    otherwise, error('Unknown quantity "%s". Use flux, bg, or survival.', quantity);
end

vars = {'longitude','latitude','energy','energy_label', zvar, 'ena_intensity_stat_uncert'};
raw  = cdfread(cdf_file, 'Variables', vars, 'CombineRecords', true);

lon_c        = double(raw{1}(:)).';     % pixel-center longitude (deg)
lat_c        = double(raw{2}(:)).';     % pixel-center latitude  (deg)
energy_kev   = double(raw{3}(:));       % energy passband centers (keV)
energy_label = raw{4};

Zraw = double(raw{5});   % selected quantity, all epochs/energies/pixels
Eraw = double(raw{6});   % ena_intensity_stat_uncert, same shape as Zraw

% Drop/select the record (epoch) dimension. cdfread with CombineRecords
% puts records along dim 1. If there's only one record, that dim may
% already be singleton/absent, so index only when it's actually present.
n_records = size(Zraw, 1);
if ndims(Zraw) == 4 && n_records > 1
    Zraw = squeeze(Zraw(epoch_idx, :, :, :));
    Eraw = squeeze(Eraw(epoch_idx, :, :, :));
else
    Zraw = squeeze(Zraw);
    Eraw = squeeze(Eraw);
end

% Auto-detect [energy, longitude, latitude] order and permute to that order
n_e = numel(energy_kev); n_o = numel(lon_c); n_a = numel(lat_c);
sz  = size(Zraw);
dim_e = find(sz == n_e, 1);
dim_o = find(sz == n_o, 1);
dim_a = find(sz == n_a, 1);
if isempty(dim_e) || isempty(dim_o) || isempty(dim_a)
    error(['Could not match %s dimensions [%s] to numel(energy)=%d, ' ...
           'numel(longitude)=%d, numel(latitude)=%d. Run size(Zraw) and ' ...
           'cdfinfo(cdf_file) to inspect manually.'], ...
           zvar, num2str(sz), n_e, n_o, n_a);
end
Zraw = permute(Zraw, [dim_e dim_o dim_a]);   % -> [energy x lon x lat]
Eraw = permute(Eraw, [dim_e dim_o dim_a]);

% Select the requested energy/ESA step, giving a [lon x lat] slice
Zsel = squeeze(Zraw(esa_step, :, :));
Esel = squeeze(Eraw(esa_step, :, :));

% Rest of the pipeline below expects Z, E as [lat x lon] (rows=lat, cols=lon)
Z = Zsel.';
E = Esel.';

% No explicit FILLVAL masking needed: fill values in this product are
% large negative numbers, and the existing Zmask logic further below
% (~isnan(Z) & Z>0) already excludes them along with true zeros/NaNs.

try
    if iscell(energy_label)
        lbl = strtrim(char(energy_label{esa_step}));
    else
        lbl = strtrim(energy_label(esa_step, :));
    end
catch
    lbl = '(unavailable)';
end
fprintf('Loaded %s\n', cdf_name);
fprintf('  quantity = %s, esa_step = %d (energy = %.4f keV, label = %s)\n', ...
    quantity, esa_step, energy_kev(esa_step), lbl);

% surfm's 'flat' shading treats grid vertex (lat(i),lon(j)) as a pixel's
% SOUTH-WEST corner and fills the cell to its north-east. The CDF gives
% pixel CENTERS, so we shift each by half a bin width (south, west) to
% turn centers into SW-corner edges. (With flip_view negating longitude,
% "west before the flip" becomes "right on screen after the flip" --
% i.e. this is the half-pixel-down / half-pixel-right correction.)
dlon = min(diff(sort(unique(lon_c))));   % longitude bin width (deg)
dlat = min(diff(sort(unique(lat_c))));   % latitude  bin width (deg)
fprintf('  bin width: dlon = %.3f deg, dlat = %.3f deg\n', dlon, dlat);

lat = lat_c - dlat/2;
lon = lon_c + dlon/2;

fprintf('  grid: %d longitude bins x %d latitude bins\n', numel(lon), numel(lat));

%==================== LONGITUDE SETUP ====================%
% Convert to [-180,180)
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

set(gca, 'ColorScale', 'linear')
colormap("turbo")

% title(sprintf('%s | E = %.3f keV (step %d)', quantity, energy_kev(esa_step), esa_step), 'FontSize', 18);
title(sprintf('%s | E = %.3f keV (step %d)', quantity, energy_kev(esa_step), esa_step), 'FontSize', 18);

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
cb.Label.String = sprintf( ...
    'Intensity (counts cm^{-2} s^{-1} sr^{-1} keV^{-1}) [ECLIPJ2000] at E = %.3f keV', ...
    energy_kev(esa_step));
else
cb.Label.String = sprintf( ...
    'Rate (counts s^{-1}) [ECLIPJ2000] at E = %.3f keV', ...
    energy_kev(esa_step));
end
cb.FontSize = 14;

if ~rate
switch esa_step
    case 1
        clim([0 2e5])
    case 2
        clim([0 1e5])
    case 3
        clim([0 3e4])
    case 4
        clim([0 6e3])
    case 5
        clim([0 1500])
    case 6
        clim([0 600])
    case 7
        clim([0 300])
end
else
    
end


%% LOS RAM PRESSURE (per pixel)
%==================== LOS RAM PRESSURE ====================%
% Implements:
%   P_ram . LOS = (2*pi*m^2*u_R^2/n_H) *
%       Integral_Emin^Emax (dE_o/E_o) j_o,ENA(E_o) (|v_o|+u_R)^2 / [sigma(E_p)*|v_o|]
%
% v_o = sqrt(2*E_o/m)                observer-frame speed
% v_p = |v_o| + u_R                  plasma-frame speed (radially inward assumption)
% E_p = m*v_p^2/2                    plasma-frame energy (keV, used only inside sigma())
% sigma(E) = Lindsay & Stebbings (2005) H+/H charge-exchange cross section,
%            valid ~0.005-600 keV:
%   sigma(E)[cm^2] = 1e-16*(1-exp(-67.3/E))^4.5*(4.15-0.531*ln(E))^2,  E in keV
%
% Implemented as a Riemann sum over the 7 discrete, non-overlapping ESA
% energy channels rather than a continuous integral: channel j
% contributes integrand(E_j) * (dE_j/E_j), where dE_j is that channel's
% full energy width (energy_delta_plus + energy_delta_minus from the
% CDF). Summing all 7 channels approximates the full-range integral;
% each channel alone gives its individual per-step contribution.
%
% ---- UNITS (worked out explicitly, in CGS: g, cm, s, erg) ----
% Dimensional analysis of the RHS: [m^2][u_R^2]/[n_H] * [j][kinematic
% bracket], with steradians dimensionless (sr=1, not a real physical
% unit) gives, after full cancellation:
%   g^2 * (cm/s)^2 / cm^-3  *  (cm^-2 s^-1 erg^-1) * (cm^-1 s^-1)
%     = g^2 cm^2 s^-4 erg^-1  =  g/s^2  =  erg/cm^2
% i.e. energy/area = pressure*length, exactly the "P_ram . LOS" the
% formula names -- PROVIDED every quantity is in truly consistent CGS,
% which requires j in cm^-2 s^-1 sr^-1 erg^-1 (per erg), NOT per keV.
%
% *** TWO BUGS FOUND AND FIXED (both confirmed against literature/an
% actual Galli et al. 2024 ApJ 971:2 comparison, whose Table 1-6 and
% Figure 9 give real ΔP*l values of order 1-500 pdyne*AU/cm^2):
%  1) j must be converted from the CDF's native per-keV to per-erg
%     (divide by keV2erg) for the CGS cancellation above to hold.
%  2) The Lindsay & Stebbings sigma(E) formula's E argument is itself
%     in keV, NOT eV as I originally (incorrectly) assumed -- verified
%     by checking sigma(1 keV) ~ 1.7e-15 cm^2 against the ~2e-15 cm^2
%     commonly quoted in the literature for H+/H at 1 keV. Passing E in
%     eV instead crushes sigma by ~1e7 (since it sits inside a
%     1/E exponential), which was the dominant source of the earlier
%     ~1e9-1e11 overestimate the person flagged after Figure-9/Table
%     comparison. ***
%
% The raw CGS result (erg/cm^2) is then converted to pdyne*AU/cm^2
% (picodyne * astronomical-unit / cm^2 -- a natural choice since
% dyne/cm^2 * AU = dyne*AU/cm^2 is dimensionally pressure*length, and
% "pico" keeps the numbers a convenient size) via:
%   1 erg/cm^2 = 1 dyne/cm = (1 cm)/(1 AU) dyne*AU/cm^2
%              = 1e12 * (1 cm)/(1 AU) pdyne*AU/cm^2

% ---- free/assumed model parameters: set these for your analysis ----
u_R_kms = 100;     % radial/bulk plasma speed, km/s (Galli et al. 2024's
                    % primary "dynamic" case; they also quote a
                    % stationary case at u_R=0, and 40-200 km/s as the
                    % plausible heliosheath range per Zirnstein et al. 2021)
n_H_cm3 = 0.10;     % assumed neutral H density, cm^-3

% ---- physical constants (CGS) ----
m_H_g   = 1.6726e-24;      % hydrogen atom mass, g
eV2erg  = 1.602177e-12;    % 1 eV in erg
keV2erg = 1e3*eV2erg;      % 1 keV in erg
kms2cms = 1e5;             % km/s -> cm/s
AU_cm   = 1.495978707e13;  % 1 AU in cm

erg_to_pdyneAU = 1e12 * (1/AU_cm);   % erg/cm^2 -> pdyne*AU/cm^2, ~=0.066846

u_R = u_R_kms * kms2cms;   % cm/s

sigma_LS05 = @(E_kev) 1e-16 .* (1 - exp(-67.3./E_kev)).^4.5 .* ...
                     (4.15 - 0.531.*log(E_kev)).^2;   % cm^2, E_kev in keV (NOT eV --
                     % verified against literature: sigma(1 keV) ~ 1.7e-15 cm^2,
                     % matching the ~2e-15 cm^2 commonly quoted for H+/H at 1 keV;
                     % passing E in eV here instead silently crushes sigma by ~1e7,
                     % which was the dominant source of the previous huge overestimate)



% ---- per-channel energy bin widths (keV), read from the CDF ----
% NOTE: cdfread's output for a single-variable request isn't always
% wrapped in a cell (this varies by MATLAB/CDF-library version), so
% unwrap defensively rather than assuming edp{1} always works.
edp_raw = cdfread(cdf_file, 'Variables', {'energy_delta_plus'},  'CombineRecords', true);
edm_raw = cdfread(cdf_file, 'Variables', {'energy_delta_minus'}, 'CombineRecords', true);
if iscell(edp_raw), edp_raw = edp_raw{1}; end
if iscell(edm_raw), edm_raw = edm_raw{1}; end
dE_kev = double(edp_raw(:)) + double(edm_raw(:));   % full width per ESA step, keV

% ---- per-channel kinematics (same for every pixel) ----
E_o_erg = energy_kev(:) * keV2erg;           % [7x1]
v_o     = sqrt(2*E_o_erg ./ m_H_g);          % cm/s
v_p     = v_o + u_R;                         % cm/s
E_p_kev = (0.5*m_H_g.*v_p.^2) ./ keV2erg;    % keV -- input to sigma_LS05
sigma_p = sigma_LS05(E_p_kev);               % cm^2

kin_factor = (v_p.^2) ./ (sigma_p .* v_o);   % [7x1] = (|v_o|+u_R)^2/(sigma(E_p)*|v_o|)
prefactor  = 2*pi * m_H_g^2 * u_R^2 / n_H_cm3;

% ---- per-pixel, per-channel contribution ----
% Zraw ([energy x lon x lat], cm^-2 s^-1 sr^-1 keV^-1) was built in the
% READ DATA section above and is untouched by the esa_step selection.
% Convert to per-erg here -- this is the unit fix described above.
j_full = Zraw / keV2erg;   % now cm^-2 s^-1 sr^-1 erg^-1, but still in the CDF's
                            % ORIGINAL (unsorted, unflipped) longitude order
j_full = j_full(:, isrt, :);   % reorder lon axis with the SAME isrt used above for
                                % Z (LONGITUDE SETUP section) so this matches the
                                % final sorted/flip_view'd 'lon' vector used for
                                % plotting -- without this, the pressure maps were
                                % misaligned/mirrored relative to the flux map
j_full(~(j_full > 0)) = NaN;   % drop fills/non-detections, same convention as Zmask

dP_step = nan(size(j_full));   % [7 x nlon x nlat], units: erg/cm^2 per channel
for jstep = 1:numel(energy_kev)
    dP_step(jstep,:,:) = prefactor * kin_factor(jstep) * ...
        (dE_kev(jstep)/energy_kev(jstep)) * j_full(jstep,:,:);
end

P_ram_step_erg = permute(dP_step, [3 2 1]);      % [lat x lon x 7], erg/cm^2 -- per-ESA-step maps
P_ram_full_erg = squeeze(nansum(dP_step, 1)).';  % [lat x lon], erg/cm^2      -- summed over all 7

P_ram_step = P_ram_step_erg * erg_to_pdyneAU;    % [lat x lon x 7], pdyne*AU/cm^2
P_ram_full = P_ram_full_erg * erg_to_pdyneAU;    % [lat x lon],      pdyne*AU/cm^2

fprintf('LOS ram pressure: u_R = %g km/s, n_H = %g cm^-3, output units: pdyne*AU/cm^2\n', u_R_kms, n_H_cm3);

%---- plot each ESA step, then the full-range sum, in the same Mollweide/ribbon style ----%
for jstep = 1:numel(energy_kev)
    plot_ena_mollweide(P_ram_step(:,:,jstep), lat, lon, flip_view, lat_center, lon_center, ...
        sprintf('LOS Ram Pressure | E = %.3f keV (step %d)', energy_kev(jstep), jstep), ...
        'P_{ram}\cdotLOS (pdyne\cdotAU/cm^2)', 20+jstep);
end
plot_ena_mollweide(P_ram_full, lat, lon, flip_view, lat_center, lon_center, ...
    'LOS Ram Pressure | Full energy range (sum of 7 steps)', ...
    'P_{ram}\cdotLOS (pdyne\cdotAU/cm^2)', 30);


%% RIBBON vs REST SUMS (per ESA step)
%==================== RIBBON vs REST SUMS ====================%
% Classify every pixel as inside vs outside the IBEX ribbon annulus
% (same Funsten et al. 2010 weighted-center geometry as the overlay
% drawn below: center + inner/outer radii) and compare flux and LOS
% ram pressure inside vs outside, per ESA step and for the full range.
%
% NOTE ON "SUM": a plain pixel-count sum will almost always be
% dominated by the "rest of sky" bucket simply because it has far more
% pixels than the thin ribbon annulus -- that's a pixel-count effect,
% not necessarily a physical one. So alongside the raw sum, this also
% reports the per-pixel MEAN (fairer for "is the ribbon enhanced"), the
% pixel counts themselves, and a solid-angle-weighted sum using the
% CDF's own 'solid_angle' variable (exact sr per pixel) as the most
% physically meaningful "total in this region" estimate.

% ---- ribbon geometry (matches the overlay in IBEX RIBBON GUIDE below) ----
rib_lon0  = 218.33;   % deg, Funsten et al. 2010 weighted-center version
rib_lat0  =  40.38;   % deg
rib_width = 36;       % deg
rib_rmid  = 74.81;    % deg
rib_rin   = rib_rmid - rib_width/2;
rib_rout  = rib_rmid + rib_width/2;

rib_lon0_plot = rib_lon0;
if flip_view
    rib_lon0_plot = -rib_lon0_plot;
end
rib_lon0_plot = mod(rib_lon0_plot + 180, 360) - 180;   % match 'lon' convention

% ---- angular distance of every pixel from the ribbon center ----
[LonG, LatG] = meshgrid(lon, lat);   % [nlat x nlon], same orientation as Z/P_ram_step
rho = acosd( sind(LatG).*sind(rib_lat0) + ...
             cosd(LatG).*cosd(rib_lat0).*cosd(LonG - rib_lon0_plot) );

ribbon_mask = (rho >= rib_rin) & (rho <= rib_rout);
rest_mask   = ~ribbon_mask;

% ---- pull per-pixel solid angle (sr) from the CDF for weighted sums ----
try
    sa_raw = cdfread(cdf_file, 'Variables', {'solid_angle'}, 'CombineRecords', true);
    if iscell(sa_raw), sa_raw = sa_raw{1}; end
    sa_raw = double(squeeze(sa_raw));               % native [lon x lat] or [lat x lon], unknown order yet
    % auto-orient to [lat x lon] the same way Zraw was auto-oriented
    if isequal(size(sa_raw), [numel(lon_c) numel(lat_c)])
        solid_angle_map = sa_raw.';                 % -> [lat x lon]
    else
        solid_angle_map = sa_raw;                   % already [lat x lon]
    end
    solid_angle_map = solid_angle_map(:, isrt);      % match the same lon reorder as everything else
    have_solid_angle = true;
catch
    warning('Could not read solid_angle from CDF; skipping solid-angle-weighted sums.');
    have_solid_angle = false;
end

% ---- per-step flux, native units, correctly lon-ordered ----
Flux_all_steps = permute(Zraw(:, isrt, :), [3 2 1]);   % [lat x lon x 7], cm^-2 s^-1 sr^-1 keV^-1
Flux_all_steps(~(Flux_all_steps > 0)) = NaN;
Flux_full = nansum(Flux_all_steps, 3);                 % [lat x lon], summed over all 7 steps

n_rib  = sum(ribbon_mask(:));
n_rest = sum(rest_mask(:));

fprintf('\n==== Ribbon vs rest (ribbon = %.1f-%.1f deg from [lat0=%.2f, lon0=%.2f]) ====\n', ...
    rib_rin, rib_rout, rib_lat0, rib_lon0_plot);
fprintf('Pixel counts: ribbon = %d, rest = %d\n', n_rib, n_rest);
fprintf('%-6s | %11s %11s %8s | %11s %11s %8s\n', ...
    'Step', 'Flux sum', 'Flux mean', 'F ratio', 'Pram sum', 'Pram mean', 'P ratio');

for jstep = 1:numel(energy_kev)
    Fmap = Flux_all_steps(:,:,jstep);
    Pmap = P_ram_step(:,:,jstep);

    f_rib_sum = nansum(Fmap(ribbon_mask));  f_rib_mean = nanmean(Fmap(ribbon_mask));
    f_rest_sum = nansum(Fmap(rest_mask));   f_rest_mean = nanmean(Fmap(rest_mask));
    p_rib_sum = nansum(Pmap(ribbon_mask));  p_rib_mean = nanmean(Pmap(ribbon_mask));
    p_rest_sum = nansum(Pmap(rest_mask));   p_rest_mean = nanmean(Pmap(rest_mask));

    fprintf('%-6d | %11.4g %11.4g %8.3f | %11.4g %11.4g %8.3f\n', ...
        jstep, f_rib_sum, f_rib_mean, f_rib_mean/f_rest_mean, ...
        p_rib_sum, p_rib_mean, p_rib_mean/p_rest_mean);
end

f_rib_sum = nansum(Flux_full(ribbon_mask));  f_rib_mean = nanmean(Flux_full(ribbon_mask));
f_rest_sum = nansum(Flux_full(rest_mask));   f_rest_mean = nanmean(Flux_full(rest_mask));
p_rib_sum = nansum(P_ram_full(ribbon_mask)); p_rib_mean = nanmean(P_ram_full(ribbon_mask));
p_rest_sum = nansum(P_ram_full(rest_mask));  p_rest_mean = nanmean(P_ram_full(rest_mask));
fprintf('%-6s | %11.4g %11.4g %8.3f | %11.4g %11.4g %8.3f\n', ...
    'Full', f_rib_sum, f_rib_mean, f_rib_mean/f_rest_mean, ...
    p_rib_sum, p_rib_mean, p_rib_mean/p_rest_mean);

if have_solid_angle
    fprintf('\n-- solid-angle-weighted (sum of value*sr, more physically meaningful than a raw pixel sum) --\n');
    fprintf('%-6s | %11s | %11s\n', 'Step', 'Flux*sr', 'Pram*sr');
    for jstep = 1:numel(energy_kev)
        Fmap = Flux_all_steps(:,:,jstep);
        Pmap = P_ram_step(:,:,jstep);
        f_rib_sa = nansum(Fmap(ribbon_mask).*solid_angle_map(ribbon_mask));
        p_rib_sa = nansum(Pmap(ribbon_mask).*solid_angle_map(ribbon_mask));
        fprintf('%-6d | %11.4g | %11.4g   (ribbon only, vs rest below)\n', jstep, f_rib_sa, p_rib_sa);
    end
    f_rest_sa_full = nansum(Flux_full(rest_mask).*solid_angle_map(rest_mask));
    p_rest_sa_full = nansum(P_ram_full(rest_mask).*solid_angle_map(rest_mask));
    fprintf('Full range: ribbon Flux*sr = %.4g, rest Flux*sr = %.4g\n', ...
        nansum(Flux_full(ribbon_mask).*solid_angle_map(ribbon_mask)), f_rest_sa_full);
    fprintf('Full range: ribbon Pram*sr = %.4g, rest Pram*sr = %.4g\n', ...
        nansum(P_ram_full(ribbon_mask).*solid_angle_map(ribbon_mask)), p_rest_sa_full);
end

%---- quick visual: per-step ribbon vs rest mean, flux and pressure ----%
rib_flux_mean_step  = nan(1, numel(energy_kev));
rest_flux_mean_step = nan(1, numel(energy_kev));
rib_pram_mean_step  = nan(1, numel(energy_kev));
rest_pram_mean_step = nan(1, numel(energy_kev));
for jstep = 1:numel(energy_kev)
    Fmap = Flux_all_steps(:,:,jstep);
    Pmap = P_ram_step(:,:,jstep);
    rib_flux_mean_step(jstep)  = nanmean(Fmap(ribbon_mask));
    rest_flux_mean_step(jstep) = nanmean(Fmap(rest_mask));
    rib_pram_mean_step(jstep)  = nanmean(Pmap(ribbon_mask));
    rest_pram_mean_step(jstep) = nanmean(Pmap(rest_mask));
end

figure(40); set(gcf,'Color','w');
subplot(1,2,1)
bar(categorical(string(round(energy_kev,3))), [rib_flux_mean_step(:) rest_flux_mean_step(:)])
legend('Ribbon (mean)','Rest of sky (mean)','Location','best')
ylabel('Intensity (cm^{-2} s^{-1} sr^{-1} keV^{-1})')
xlabel('Energy (keV)')
title('Flux: ribbon vs rest, per ESA step')
set(gca,'YScale','log'); grid on

subplot(1,2,2)
bar(categorical(string(round(energy_kev,3))), [rib_pram_mean_step(:) rest_pram_mean_step(:)])
legend('Ribbon (mean)','Rest of sky (mean)','Location','best')
ylabel('P_{ram}\cdotLOS (pdyne\cdotAU/cm^2)')
xlabel('Energy (keV)')
title('LOS Ram Pressure: ribbon vs rest, per ESA step')
set(gca,'YScale','log'); grid on

figure(1)


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
set(gca,'YScale','linear')
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
set(gca,'YScale','linear')
xlabel('Angular distance across ribbon (deg)')
ylabel('Average Intensity')
title('Stacked Ribbon Profile (Average)')
grid on

figure(4)
errorbar(s_cut, ave_profile, ave_profile_err, 'k-', 'LineWidth', 1.5)
set(gca,'YScale','linear')
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
set(gca, 'ColorScale', 'linear')

figure(1)

if save
    T = table(s_cut(:), ave_profile(:), ave_profile_err(:),...
        'VariableNames', {'s_cut_deg','average_intensity','ave_profile_err'});
    writetable(T, [fld quantity '_ave_esa' num2str(esa_step) '.csv']);
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
set(gca,'YScale','linear')
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
set(gca,'ColorScale','linear')
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
    'Intensity (counts cm^{-2} s^{-1} sr^{-1} keV^{-1}) [ECLIPJ2000] at E = %.3f keV', ...
    energy_kev(esa_step));
cb.FontSize = 14;




%% Functions

function plot_ena_mollweide(Zmap, lat, lon, flip_view, lat_center, lon_center, ...
    title_str, cbar_label, fig_num)
% Mollweide map + IBEX ribbon guide/markers, factored out of the main
% script so it can be reused for the LOS ram pressure maps (same
% projection, wrap-padding, and overlay style as the primary flux plot;
% does NOT include the cross-ribbon-cut/profile analysis). No clim() is
% applied here since the pressure ranges are model-dependent -- add one
% manually if you want to fix the color scale across figures.
%
% Zmap, lat, lon must already be on the same [lat x lon] grid / edge
% convention used elsewhere in this script (i.e. lat/lon already shifted
% to SW-corner edges, see the half-pixel correction earlier).

lon_ext = [lon, lon(1) + 360];
Z_ext   = [Zmap, Zmap(:,1)];
lat_ext = [lat, lat(end) + (lat(2)-lat(1))];
Z_ext   = [Z_ext; Z_ext(end,:)];

[Lon, Lat] = meshgrid(lon_ext, lat_ext);
Zmask = ~isnan(Z_ext) & (Z_ext ~= 0);
Z_ext(~Zmask) = NaN;

figure(fig_num);
set(gcf, 'Color', [1 1 1]);
ax = axesm('mollweid', ...
    'Origin',[lat_center lon_center 0], ...
    'MapLatLimit', [-90 90], ...
    'MapLonLimit', [-180 180], ...
    'Frame', 'on', ...
    'Grid', 'on', ...
    'MeridianLabel', 'off', ...
    'ParallelLabel', 'off', ...
    'MLineLocation', 30, ...
    'PLineLocation', 30, ...
    'GLineWidth', 1.0, ...
    'Gcolor', [1 0.6 0], ...
    'FEdgeColor', [0 0 0.5], ...
    'FontSize', 12);
axis off
tightmap
set(ax, 'Color', 'none')
setm(ax, 'FFaceColor', [0 0 0])

h = surfm(Lat, Lon, Z_ext);
set(h, 'EdgeColor', 'none', 'FaceColor', 'flat', ...
    'FaceAlpha', 'flat', 'AlphaData', double(Zmask), 'AlphaDataMapping', 'none');
set(gca, 'ColorScale', 'linear')
colormap("turbo")
title(title_str, 'FontSize', 16);

% Manual meridian labels
lablat = 0;
meridians = -180:30:180;
for k = 1:numel(meridians)
    m = meridians(k);
    val = mod(-m + 180, 360) - 180;
    if abs(val) == 180
        txt = '180°';
    elseif val == 0
        txt = '0°';
    else
        txt = sprintf('%d°', val);
    end
    textm(lablat, m, txt, 'Color', [1 1 1], 'FontSize', 12, ...
        'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end

cb = colorbar('southoutside');
cb.Label.String = cbar_label;
cb.FontSize = 12;

% IBEX ribbon guide (Funsten et al. 2010 weighted-center version)
lon0 = 218.33; lat0 = 40.38; width = 36;
r_mid = 74.81; r_in = r_mid - width/2; r_out = r_mid + width/2;
az = linspace(0,360,721);
[lat_mid, lon_mid] = smallcircle_deg(lat0, lon0, r_mid, az);
[lat_in , lon_in ] = smallcircle_deg(lat0, lon0, r_in , az);
[lat_out, lon_out] = smallcircle_deg(lat0, lon0, r_out, az);

lon_v1 = mod(255.0 + 180, 360) - 180; lat_v1 =  35.0;
lon_v2 = mod(290.4 + 180, 360) - 180; lat_v2 = -32.2;
lon_nose = mod(259.200622 + 180, 360) - 180; lat_nose =  5.116296;
lon_tail = mod( 79.200622 + 180, 360) - 180; lat_tail = -5.116296;

if flip_view
    lon_mid = -lon_mid; lon_in = -lon_in; lon_out = -lon_out;
    lon_v1 = -lon_v1; lon_v2 = -lon_v2;
    lon_nose = -lon_nose; lon_tail = -lon_tail; lon0 = -lon0;
end
lon_mid = mod(lon_mid + 180, 360) - 180;
lon_in  = mod(lon_in  + 180, 360) - 180;
lon_out = mod(lon_out + 180, 360) - 180;
lon_v1  = mod(lon_v1  + 180, 360) - 180;
lon_v2  = mod(lon_v2  + 180, 360) - 180;
lon_nose= mod(lon_nose+ 180, 360) - 180;
lon_tail= mod(lon_tail+ 180, 360) - 180;
lon0    = mod(lon0    + 180, 360) - 180;

plotm(lat_mid, lon_mid, 'w-',  'LineWidth', 2)
plotm(lat_in , lon_in , 'w--', 'LineWidth', 1.5)
plotm(lat_out, lon_out, 'w--', 'LineWidth', 1.5)

plotm(lat0, lon0, 'wo', 'MarkerFaceColor','w', 'MarkerSize',6)
textm(lat0, lon0, 'RC', 'Color','w','FontSize',12,'FontWeight','bold', ...
    'HorizontalAlignment','right','VerticalAlignment','bottom')

plotm(lat_v1, lon_v1, 'wo', 'MarkerFaceColor','w', 'MarkerSize',6)
textm(lat_v1, lon_v1, '  V1', 'Color','w','FontSize',12,'FontWeight','bold',...
    'HorizontalAlignment','right', 'VerticalAlignment','bottom')

plotm(lat_v2, lon_v2, 'wo', 'MarkerFaceColor','w', 'MarkerSize',6)
textm(lat_v2, lon_v2, '  V2', 'Color','w','FontSize',12,'FontWeight','bold',...
    'HorizontalAlignment','right', 'VerticalAlignment','top')

plotm(lat_nose, lon_nose, 'wo', 'MarkerFaceColor','w', 'MarkerSize',6)
textm(lat_nose, lon_nose-3, '  Nose', 'Color','w','FontSize',12,'FontWeight','bold',...
    'HorizontalAlignment','left', 'VerticalAlignment','bottom')

plotm(lat_tail, lon_tail, 'wo', 'MarkerFaceColor','w', 'MarkerSize',6)
textm(lat_tail, lon_tail, '  Tail', 'Color','w','FontSize',12,'FontWeight','bold',...
    'HorizontalAlignment','right', 'VerticalAlignment','top')

end

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