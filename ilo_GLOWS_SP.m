%% IMAP GLOWS L3e
% 1) Scalar variables: plot versus time, markers only, no repoint text labels
% 2) 2D variables: still plot each 2D variable for each CDF file
% 3) 3D+ variables: plot first slice only

close all; clearvars; clc; format long; warning('off','all')

rootdir = ['/Users/mitchellshen/Library/CloudStorage/Dropbox-Personal/' ...
    'IMAP_ILO/SDC/20260526_SP_corr_maps/input_glows/l3e'];

outdir_scalar = fullfile(rootdir, 'fig_scalar_time_series');
outdir_2d     = fullfile(rootdir, 'fig_2d_vars');

if ~exist(outdir_scalar, 'dir'); mkdir(outdir_scalar); end
if ~exist(outdir_2d, 'dir'); mkdir(outdir_2d); end

files = dir(fullfile(rootdir, '**', '*.cdf'));
nFiles = numel(files);
fprintf('Found %d CDF files\n', nFiles);

dates   = NaT(nFiles,1);
repoint = nan(nFiles,1);

S = struct();

%% Loop over files

for i = 1:nFiles

    fpath = fullfile(files(i).folder, files(i).name);
    fprintf('\n[%d/%d] %s\n', i, nFiles, files(i).name);

    tok = regexp(files(i).name, '_(\d{8})-repoint(\d+)_', 'tokens', 'once');

    if ~isempty(tok)
        date_str = tok{1};
        repoint_str = tok{2};

        dates(i)   = datetime(date_str, 'InputFormat','yyyyMMdd');
        repoint(i) = str2double(repoint_str);
    else
        date_str = sprintf('file%03d', i);
        repoint_str = sprintf('%05d', i);

        dates(i)   = NaT;
        repoint(i) = i;
    end

    info = cdfinfo(fpath);
    vars = info.Variables(:,1);

    for v = 1:numel(vars)

        varname = vars{v};
        safevar = matlab.lang.makeValidName(varname);

        try
            tmp = cdfread(fpath, ...
                'Variables', {varname}, ...
                'ConvertEpochToDatenum', true);

            data = tmp{1};

            if ~isnumeric(data)
                continue
            end

            data = double(squeeze(data));

            if isempty(data) || all(isnan(data(:)))
                continue
            end

            sz = size(data);
            fprintf('  %s size = %s\n', varname, mat2str(sz));

            %% Scalar variables: collect for later time-series plotting

            if isscalar(data)

                if ~isfield(S, safevar)
                    S.(safevar).original_name = varname;
                    S.(safevar).value = nan(nFiles,1);
                end

                S.(safevar).value(i) = data;

            %% 2D variables: plot immediately

            elseif ismatrix(data)

                fig = figure('Color','w', ...
                    'Position',[100 100 900 650], ...
                    'Visible','off');

                imagesc(data);
                axis xy;
                colorbar;

                xlabel('Column index');
                ylabel('Row index');

                title(sprintf('%s | %s | repoint %s | size=%s', ...
                    varname, date_str, repoint_str, mat2str(sz)), ...
                    'Interpreter','none');

                set(gca, 'FontSize', 14, 'LineWidth', 1.2);

                fname_out = sprintf('%s_repoint%s_%s_2D.png', ...
                    date_str, repoint_str, safevar);

                exportgraphics(fig, fullfile(outdir_2d, fname_out), ...
                    'Resolution', 200);

                close(fig);

            %% 3D or higher: plot first slice

            else

                fig = figure('Color','w', ...
                    'Position',[100 100 900 650], ...
                    'Visible','off');

                data2 = squeeze(data(:,:,1));

                imagesc(data2);
                axis xy;
                colorbar;

                xlabel('Column index');
                ylabel('Row index');

                title(sprintf('%s first slice | %s | repoint %s | size=%s', ...
                    varname, date_str, repoint_str, mat2str(sz)), ...
                    'Interpreter','none');

                set(gca, 'FontSize', 14, 'LineWidth', 1.2);

                fname_out = sprintf('%s_repoint%s_%s_firstSlice.png', ...
                    date_str, repoint_str, safevar);

                exportgraphics(fig, fullfile(outdir_2d, fname_out), ...
                    'Resolution', 200);

                close(fig);
            end

        catch ME
            fprintf('  skipped %s: %s\n', varname, ME.message);
        end
    end
end

%% Sort scalar variables by date / repoint

valid_date = ~isnat(dates);
sort_key_date = datenum(dates);
sort_key_date(~valid_date) = inf;

[~, idx] = sortrows([sort_key_date(:), repoint(:)]);

dates   = dates(idx);
repoint = repoint(idx);

fields = fieldnames(S);

for k = 1:numel(fields)
    S.(fields{k}).value = S.(fields{k}).value(idx);
end

%% Plot scalar variables versus time
% Markers only. No repoint labels.

for k = 1:numel(fields)

    safevar = fields{k};
    varname = S.(safevar).original_name;
    y = S.(safevar).value;

    good = ~isnan(y) & ~isnat(dates);

    if nnz(good) < 1
        continue
    end

    fig = figure('Color','w', ...
        'Position',[100 100 1000 650], ...
        'Visible','off');

    plot(dates(good), y(good), 'o', ...
        'LineWidth', 1.5, ...
        'MarkerSize', 6);

    grid on;
    xlabel('Date');
    ylabel(varname, 'Interpreter','none');

    title(sprintf('%s versus time', varname), ...
        'Interpreter','none');

    xtickformat('yyyy-MM-dd');
    set(gca, 'FontSize', 14, 'LineWidth', 1.2);

    fname_out = sprintf('time_series_%s.png', safevar);

    exportgraphics(fig, fullfile(outdir_scalar, fname_out), ...
        'Resolution', 200);

    close(fig);
end

%% Optional combined spacecraft figure

wanted = {'spacecraft_longitude', ...
          'spacecraft_latitude', ...
          'spacecraft_radius', ...
          'spacecraft_velocity_x', ...
          'spacecraft_velocity_y', ...
          'spacecraft_velocity_z', ...
          'spin_axis_longitude', ...
          'spin_axis_latitude', ...
          'elongation', ...
          'glows_flags'};

fig = figure('Color','w', ...
    'Position',[100 100 1300 850], ...
    'Visible','off');

tiledlayout(2,5, 'TileSpacing','compact', 'Padding','compact');

for q = 1:numel(wanted)

    nexttile;

    safevar = matlab.lang.makeValidName(wanted{q});

    if isfield(S, safevar)

        y = S.(safevar).value;
        good = ~isnan(y) & ~isnat(dates);

        plot(dates(good), y(good), 'o', ...
            'LineWidth', 1.2, ...
            'MarkerSize', 4);

        title(wanted{q}, 'Interpreter','none');
        grid on;
        xtickformat('yyyy-MM');
        set(gca, 'FontSize', 10, 'LineWidth', 1.0);

    else
        text(0.5, 0.5, 'not found', ...
            'HorizontalAlignment','center');
        axis off;
        title(wanted{q}, 'Interpreter','none');
    end
end

sgtitle('Scalar variables versus time');

exportgraphics(fig, fullfile(outdir_scalar, ...
    'combined_scalar_time_series_markers_only.png'), ...
    'Resolution', 200);

close(fig);

fprintf('\nDone.\n');
fprintf('Scalar time-series figures saved to:\n%s\n', outdir_scalar);
fprintf('2D / first-slice figures saved to:\n%s\n', outdir_2d);