function [matched_params, analysis_data] = multiband_spectrograph_sweep(band_selection, varargin)
% MULTIBAND_SPECTROGRAPH_SWEEP - Sweep Y/J/H spectrograph parameters
%
% Calculates how focal lengths, beam diameter, and projected grating width
% vary with grating line density for selected photometric bands. A reference
% grating density is selected for the first band, and the remaining bands
% are matched to the closest camera focal length.
%
% Inputs:
%   band_selection              - {'Y','J','H'} or a subset
%   'name'                      - Design name
%   'resolving_power'           - Resolving power for each band
%   'lambda_central'            - Central wavelengths [m]
%   'slit_width'                - Entrance slit width [m]
%   'slit_image_width'          - Required slit-image widths [m]
%   'f_number_coll'             - Collimator f-number
%   'rho_range'                 - Grating density range [lines/m]
%   'reference_grating_density' - Reference density [lines/m]
%   'save_plots'                - Save figures to disk
%   'results_dir'               - Output directory for figures
%
% Outputs:
%   matched_params - Matched design parameters for the selected bands
%   analysis_data  - Complete sweep arrays and relative camera-focal spread
%
% Example:
%   [matched, data] = multiband_spectrograph_sweep({'Y','J','H'});

%% Input parsing
p = inputParser;
addRequired(p, 'band_selection', @iscell);
addParameter(p, 'name', 'MCIFU_5000_950', @(x) ischar(x) || isstring(x));
addParameter(p, 'resolving_power', [5000, 5000, 5000], @isnumeric);
addParameter(p, 'lambda_central', ...
    [1.0375e-6, 1.2475e-6, 1.635e-6], @isnumeric);
addParameter(p, 'slit_width', 7.3e-6, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'slit_image_width', ...
    [36e-6, 38e-6, 50e-6], @isnumeric);
addParameter(p, 'f_number_coll', 4.55, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'rho_range', ...
    linspace(400e3, 1200e3, 5000), @isnumeric);
addParameter(p, 'reference_grating_density', 650e3, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'save_plots', true, ...
    @(x) islogical(x) || ismember(x, [0, 1]));
addParameter(p, 'results_dir', 'Results', ...
    @(x) ischar(x) || isstring(x));

parse(p, band_selection, varargin{:});
params = p.Results;

band_map = containers.Map({'Y','J','H'}, [1,2,3]);

try
    selected_indices = cellfun(@(x) band_map(char(x)), band_selection);
catch
    error('band_selection must contain only Y, J, and/or H.');
end

bands = band_selection;

if numel(params.resolving_power) == numel(bands)
    R = params.resolving_power;
else
    R = params.resolving_power(selected_indices);
end

if numel(params.lambda_central) == numel(bands)
    lambdaC = params.lambda_central;
else
    lambdaC = params.lambda_central(selected_indices);
end

if numel(params.slit_image_width) == numel(bands)
    s2 = params.slit_image_width;
else
    s2 = params.slit_image_width(selected_indices);
end

s1 = params.slit_width;
F1 = params.f_number_coll;
rhoRange = params.rho_range;

n_bands = numel(bands);
n_points = numel(rhoRange);

fprintf('Analyzing spectrograph parameters for bands: %s\n', ...
    strjoin(bands, ', '));
fprintf('Parameter sweep range: %.0f to %.0f lines/mm (%d points)\n', ...
    min(rhoRange)*1e-3, max(rhoRange)*1e-3, n_points);

%% Core parameter calculations
% The slit-image width ratio sets the camera/collimator f-number ratio.
F2 = (s2 ./ s1) .* F1;

deltaLambda = zeros(1, n_bands);
alpha = NaN(n_bands, n_points);
f2 = NaN(n_bands, n_points);
D = NaN(n_bands, n_points);
f1 = NaN(n_bands, n_points);
W = NaN(n_bands, n_points);

for i = 1:n_bands

    deltaLambda(i) = lambdaC(i) / R(i);

    % Littrow grating equation.
    asin_arg = (rhoRange .* lambdaC(i)) / 2;
    valid = abs(asin_arg) <= 1;
    alpha(i, valid) = asin(asin_arg(valid));

    % Camera focal length from the required spectral sampling.
    f2(i, :) = ...
        (cos(alpha(i, :)) .* s2(i)) ./ ...
        (rhoRange .* deltaLambda(i));

    % In Littrow the incident and diffracted beam diameters are equal.
    D(i, :) = f2(i, :) ./ F2(i);
    f1(i, :) = F1 .* D(i, :);

    % Illuminated width projected onto the grating.
    W(i, :) = D(i, :) ./ cos(alpha(i, :));
end

%% Parameter plots
colors = lines(max(3, n_bands));

for i = 1:n_bands
    fig = create_parameter_plot(...
        bands{i}, rhoRange, f1(i,:), f2(i,:), D(i,:), W(i,:), colors);

    if params.save_plots
        save_parameter_plot(fig, params.name, bands{i}, params.results_dir);
    end
end

%% Cross-band matching
ref_band_idx = 1;

[matched_params, analysis_data] = match_cross_band_parameters(...
    bands, rhoRange, params.reference_grating_density, ...
    f1, f2, D, W, alpha, lambdaC, ref_band_idx);

display_matched_results(matched_params, bands);

analysis_data.rho_range = rhoRange;
analysis_data.f1_matrix = f1;
analysis_data.f2_matrix = f2;
analysis_data.D_matrix = D;
analysis_data.W_matrix = W;
analysis_data.alpha_matrix = alpha;
analysis_data.F2 = F2;
analysis_data.resolving_power = R;
analysis_data.lambda_central = lambdaC;

fprintf('Parameter sweep completed successfully.\n');

end

function fig = create_parameter_plot(...
    band_name, rhoRange, f1_band, f2_band, D_band, W_band, colors)

    fig = figure('Name', ...
        sprintf('Spectrograph Parameters - Band %s', band_name), ...
        'Position', [100, 100, 800, 600]);

    hold on;

    plot(rhoRange * 1e-3, f2_band * 1e3, ...
        'LineWidth', 2.5, 'Color', colors(1,:), ...
        'DisplayName', 'f_2 (camera)');

    plot(rhoRange * 1e-3, D_band * 1e3, ...
        'LineWidth', 2.5, 'Color', colors(2,:), ...
        'DisplayName', 'D (beam)');

    plot(rhoRange * 1e-3, f1_band * 1e3, ...
        'LineWidth', 2.5, 'Color', colors(3,:), ...
        'DisplayName', 'f_1 (collimator)');

    if size(colors,1) >= 4
        w_color = colors(4,:);
    else
        w_color = [0.4940, 0.1840, 0.5560];
    end

    plot(rhoRange * 1e-3, W_band * 1e3, ...
        '--', 'LineWidth', 1.8, 'Color', w_color, ...
        'DisplayName', 'W (projected grating width)');

    xlabel('Grating Density [lines/mm]');
    ylabel('Length [mm]');
    title(sprintf('Spectrograph Parameters - Band %s', band_name));

    legend('show', 'Location', 'best');
    grid on;
    set(gca, 'YScale', 'log');
end

function [matched_params, analysis_data] = match_cross_band_parameters(...
    bands, rhoRange, rho_ref_requested, f1, f2, D, W, alpha, lambdaC, ref_band_idx)

    matched_params = struct();
    n_bands = numel(bands);

    [~, idx_ref] = min(abs(rhoRange - rho_ref_requested));

    if ~isfinite(f2(ref_band_idx, idx_ref))
        valid_idx = find(isfinite(f2(ref_band_idx, :)));
        if isempty(valid_idx)
            error('No physically valid grating density found for reference band %s.', ...
                bands{ref_band_idx});
        end
        [~, k] = min(abs(rhoRange(valid_idx) - rho_ref_requested));
        idx_ref = valid_idx(k);
    end

    rho_ref = rhoRange(idx_ref);

    matched_params.ref_band = bands{ref_band_idx};
    matched_params.rho(ref_band_idx) = rho_ref;
    matched_params.f1(ref_band_idx) = f1(ref_band_idx, idx_ref);
    matched_params.f2(ref_band_idx) = f2(ref_band_idx, idx_ref);
    matched_params.D(ref_band_idx) = D(ref_band_idx, idx_ref);
    matched_params.W(ref_band_idx) = W(ref_band_idx, idx_ref);
    matched_params.alpha_deg(ref_band_idx) = ...
        rad2deg(alpha(ref_band_idx, idx_ref));

    fprintf('\n--- Cross-Band Matching ---\n');
    fprintf('Reference band %s at rho = %.0f lines/mm:\n', ...
        bands{ref_band_idx}, rho_ref*1e-3);

    fprintf('  f1 = %.2f mm, f2 = %.2f mm, D = %.2f mm, alpha = %.1f deg\n', ...
        matched_params.f1(ref_band_idx)*1e3, ...
        matched_params.f2(ref_band_idx)*1e3, ...
        matched_params.D(ref_band_idx)*1e3, ...
        matched_params.alpha_deg(ref_band_idx));

    for i = 1:n_bands

        if i == ref_band_idx
            continue;
        end

        mismatch = abs(f2(i,:) - matched_params.f2(ref_band_idx));
        mismatch(~isfinite(mismatch)) = Inf;

        [min_mismatch, idx_match] = min(mismatch);

        if ~isfinite(min_mismatch)
            error('No physically valid grating density found for band %s.', bands{i});
        end

        matched_params.rho(i) = rhoRange(idx_match);
        matched_params.f1(i) = f1(i, idx_match);
        matched_params.f2(i) = f2(i, idx_match);
        matched_params.D(i) = D(i, idx_match);
        matched_params.W(i) = W(i, idx_match);
        matched_params.alpha_deg(i) = rad2deg(alpha(i, idx_match));

        fprintf('Band %s matched to f2 = %.2f mm:\n', ...
            bands{i}, matched_params.f2(i)*1e3);

        fprintf('  rho = %.0f lines/mm, f1 = %.2f mm, D = %.2f mm, alpha = %.1f deg\n', ...
            matched_params.rho(i)*1e-3, ...
            matched_params.f1(i)*1e3, ...
            matched_params.D(i)*1e3, ...
            matched_params.alpha_deg(i));
    end

    analysis_data.relative_f2_spread = ...
        std(matched_params.f2) / mean(matched_params.f2);
end

function save_parameter_plot(fig, design_name, band_name, results_dir)

    results_dir = char(results_dir);

    if ~exist(results_dir, 'dir')
        mkdir(results_dir);
    end

    fig_name = sprintf('parameters_%s_%s', char(design_name), band_name);

    saveas(fig, fullfile(results_dir, [fig_name, '.fig']));
    saveas(fig, fullfile(results_dir, [fig_name, '.png']));
end

function display_matched_results(matched_params, bands)

    fprintf('\n=== MATCHED PARAMETERS ===\n');

    fprintf('%-6s %-12s %-10s %-10s %-10s %-10s %-12s\n', ...
        'Band', 'rho [l/mm]', 'f1 [mm]', 'f2 [mm]', ...
        'D [mm]', 'W [mm]', 'alpha [deg]');

    fprintf('%s\n', repmat('-', 76, 1));

    for i = 1:numel(bands)

        fprintf('%-6s %-12.0f %-10.1f %-10.1f %-10.1f %-10.1f %-12.1f\n', ...
            bands{i}, ...
            matched_params.rho(i) * 1e-3, ...
            matched_params.f1(i) * 1e3, ...
            matched_params.f2(i) * 1e3, ...
            matched_params.D(i) * 1e3, ...
            matched_params.W(i) * 1e3, ...
            matched_params.alpha_deg(i));
    end
end
