function [T_sorted, best_design, results] = slit_limited_spectrograph_sweep(varargin)
% SLIT_LIMITED_SPECTROGRAPH_SWEEP - Sweep slit-limited spectrograph designs
%
% Explores first-order spectrograph configurations over resolving power,
% slit width, grating density, detector format, and Littrow/off-Littrow
% geometry. Candidate designs are checked for detector coverage, camera
% focal length, and camera focal ratio.
%
% A scalar parameter set can also be used to evaluate a single design.
% If 'camera_focal_length' is empty, f2 is derived from the requested
% detector sampling. If specified, the supplied f2 is used directly.
%
% Inputs (parameter/value pairs):
%   'wavelength_range'    - [min,max] wavelength [m]
%   'sampling_min'        - Required sampling [pixels/resolution element]
%   'NA'                  - Input numerical aperture
%   'R_list'              - Resolving-power values to test
%   'slit_widths'         - Slit widths [m]
%   'grating_densities'   - Grating line densities [lines/m]
%   'detectors'           - Cell array: {pixel_um, nx, ny, name}
%   'geometries'          - {'Littrow'} or {'Littrow','OffLittrow'}
%   'off_littrow_angle'   - Incidence angle for OffLittrow [deg]
%   'camera_focal_length' - Optional fixed f2 [m]
%   'save_results'        - Save sorted table to MAT file
%   'output_file'         - MAT filename
%   'make_plot'           - Plot resolving power of highest-ranked design
%
% Outputs:
%   T_sorted     - Table of valid configurations, feasible designs first
%   best_design  - Highest-ranked feasible design as a one-row table
%   results      - Structure array with wavelength-dependent data
%
% Example:
%   [T, best, results] = slit_limited_spectrograph_sweep(...
%       'R_list', 1300:50:1500, ...
%       'grating_densities', 1511e3);

%% Input parsing
p = inputParser;

addParameter(p, 'wavelength_range', [550e-9, 700e-9], ...
    @(x) isnumeric(x) && numel(x) == 2 && all(x > 0) && x(2) > x(1));
addParameter(p, 'sampling_min', 3, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'NA', 0.11, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
addParameter(p, 'R_list', [1300, 1350, 1400, 1450, 1500], ...
    @(x) isnumeric(x) && all(x > 0));
addParameter(p, 'slit_widths', 50e-6, ...
    @(x) isnumeric(x) && all(x > 0));
addParameter(p, 'grating_densities', 1511e3, ...
    @(x) isnumeric(x) && all(x > 0));
addParameter(p, 'detectors', {5.2, 1280, 1024, 'DCC1545M'}, ...
    @(x) iscell(x) && size(x, 2) >= 4);
addParameter(p, 'geometries', {'Littrow'}, ...
    @(x) iscell(x) && all(cellfun(@(s) ischar(s) || isstring(s), x)));
addParameter(p, 'off_littrow_angle', 20, ...
    @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'camera_focal_length', [], ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
addParameter(p, 'save_results', true, ...
    @(x) islogical(x) || ismember(x, [0, 1]));
addParameter(p, 'output_file', 'sweep_results.mat', ...
    @(x) ischar(x) || isstring(x));
addParameter(p, 'make_plot', true, ...
    @(x) islogical(x) || ismember(x, [0, 1]));

parse(p, varargin{:});
params = p.Results;

lambda_min = params.wavelength_range(1);
lambda_max = params.wavelength_range(2);
lambda_c = mean(params.wavelength_range);

F1 = 1 / (2 * params.NA);

R_list = params.R_list;
s1_list = params.slit_widths;
rho_list = params.grating_densities;
detectors = params.detectors;
geometries = params.geometries;

results = struct([]);
row = 0;

%% Parameter sweep
for idet = 1:size(detectors, 1)

    pix = detectors{idet, 1} * 1e-6;
    nx = detectors{idet, 2};
    name = detectors{idet, 4};
    sensor_width = nx * pix;

    for iR = 1:numel(R_list)

        R = R_list(iR);
        delta_lambda = lambda_c / R;

        for islit = 1:numel(s1_list)

            s1 = s1_list(islit);

            for irho = 1:numel(rho_list)

                rho = rho_list(irho);

                for iga = 1:numel(geometries)

                    geom = char(geometries{iga});

                    %% Grating geometry
                    if strcmpi(geom, 'Littrow')

                        arg = rho * lambda_c / 2;
                        if abs(arg) > 1
                            continue;
                        end

                        alpha = asin(arg);
                        beta = alpha;

                    elseif strcmpi(geom, 'OffLittrow')

                        alpha = deg2rad(params.off_littrow_angle);
                        arg = rho * lambda_c - sin(alpha);

                        if abs(arg) > 1
                            continue;
                        end

                        beta = asin(arg);

                    else
                        error('Unsupported geometry: %s', geom);
                    end

                    % Verify that the full wavelength band can be diffracted.
                    arg_min = rho * lambda_min - sin(alpha);
                    arg_max = rho * lambda_max - sin(alpha);

                    if any(abs([arg_min, arg_max]) > 1)
                        continue;
                    end

                    beta_min = asin(arg_min);
                    beta_max = asin(arg_max);

                    % Beam anamorphic ratio D2/D1.
                    Ma = cos(beta) / cos(alpha);

                    %% Camera and collimator geometry
                    s2_target = params.sampling_min * pix;

                    % F2/F1 gives the slit-image magnification in the
                    % dispersion direction for the projected geometry below.
                    F2 = (s2_target / s1) * F1;

                    if isempty(params.camera_focal_length)
                        % Camera focal length required by the requested
                        % spectral sampling at the central wavelength.
                        f2 = (s2_target * cos(beta)) / (rho * delta_lambda);
                    else
                        f2 = params.camera_focal_length;
                    end

                    D2 = f2 / F2;

                    % W is the illuminated width projected onto the grating.
                    W = D2 / cos(beta);
                    D1 = W * cos(alpha);
                    f1 = F1 * D1;

                    focal_length_ratio = f2 / f1;
                    slit_magnification = ...
                        focal_length_ratio * cos(alpha) / cos(beta);

                    %% Dispersion and detector coverage
                    dlambda_dx = cos(beta) / (rho * f2);

                    spectrum_length = ...
                        abs(f2 * (tan(beta_max) - tan(beta_min)));
                    n_pixels_range = spectrum_length / pix;
                    pixels_per_res = delta_lambda / (dlambda_dx * pix);

                    %% Resolving power across wavelength
                    lambda_vec = linspace(lambda_min, lambda_max, 50);
                    R_all = (rho .* lambda_vec * F1 .* W) ./ s1;

                    %% Feasibility checks
                    fits_detector = ceil(n_pixels_range) <= nx;
                    f2_ok = (f2 >= 20e-3) && (f2 <= 1000e-3);
                    F2_ok = F2 > 1.4;

                    feasible = fits_detector && f2_ok && F2_ok;

                    %% Ranking heuristic
                    det_util = min(1, n_pixels_range / nx);
                    slit_term = s1 / max(s1_list);
                    FoM = 0.5 * det_util + 0.5 * slit_term;

                    %% Store configuration
                    row = row + 1;

                    results(row).config_id = row;
                    results(row).det_name = name;
                    results(row).pix_um = detectors{idet, 1};
                    results(row).nx = nx;
                    results(row).R = R;
                    results(row).s1_um = s1 * 1e6;
                    results(row).rho_lpm = rho * 1e-3;
                    results(row).geom = geom;
                    results(row).alpha_deg = rad2deg(alpha);
                    results(row).beta_deg = rad2deg(beta);
                    results(row).Ma = Ma;
                    results(row).f2_mm = f2 * 1e3;
                    results(row).f1_mm = f1 * 1e3;
                    results(row).D1_mm = D1 * 1e3;
                    results(row).D2_mm = D2 * 1e3;
                    results(row).F2 = F2;
                    results(row).focal_length_ratio = focal_length_ratio;
                    results(row).slit_magnification = slit_magnification;
                    results(row).dx_mm_per_nm = 1 / (dlambda_dx * 1e6);
                    results(row).pixels_per_res = pixels_per_res;
                    results(row).n_pixels_range = n_pixels_range;
                    results(row).detector_width_mm = sensor_width * 1e3;
                    results(row).fits_detector = fits_detector;
                    results(row).f2_ok = f2_ok;
                    results(row).F2_ok = F2_ok;
                    results(row).feasible = feasible;
                    results(row).FoM = FoM;
                    results(row).R_all = R_all;
                    results(row).lambda_vec = lambda_vec;

                end
            end
        end
    end
end

if isempty(results)
    error('No valid configurations were generated for the selected parameter ranges.');
end

%% Sort and optionally save results
T = struct2table(results);
T_sorted = sortrows(T, {'feasible', 'FoM'}, {'descend', 'descend'});

if params.save_results
    save(char(params.output_file), 'T_sorted');
end

%% Display table
fprintf('=== POSSIBLE DESIGNS ===\n');

fprintf('%-5s %-15s %-8s %-6s %-6s %-7s %-8s | %-5s %-6s %-6s %-5s %-7s %-6s %-5s | %-5s\n', ...
    'Rank', 'Detector', 'Pix(um)', 'Samp', 'R', 's1(um)', 'rho(l/mm)', ...
    'a(deg)', 'f1(mm)', 'f2(mm)', 'F2', 'D2(mm)', 'Mslit', '#Pix', 'FoM');

fprintf('%s\n', repmat('-', 124, 1));

n_display = min(100, height(T_sorted));

for i = 1:n_display

    feasible_char = 'Y';
    if ~T_sorted.feasible(i)
        feasible_char = 'N';
    end

    fprintf('%-5d %-15s %-8.1f %-6.1f %-6d %-7.0f %-8.0f | %-5.1f %-6.1f %-6.1f %-5.1f %-7.1f %-6.2f %-5.0f | %-5.2f %s\n', ...
        i, T_sorted.det_name{i}, T_sorted.pix_um(i), ...
        T_sorted.pixels_per_res(i), T_sorted.R(i), ...
        T_sorted.s1_um(i), T_sorted.rho_lpm(i), ...
        T_sorted.alpha_deg(i), T_sorted.f1_mm(i), ...
        T_sorted.f2_mm(i), T_sorted.F2(i), T_sorted.D2_mm(i), ...
        T_sorted.slit_magnification(i), ...
        T_sorted.n_pixels_range(i), T_sorted.FoM(i), feasible_char);
end

%% Highest-ranked feasible design
feasible_designs = T_sorted(T_sorted.feasible, :);

if isempty(feasible_designs)
    best_design = T_sorted([], :);
    fprintf('\nNo feasible designs satisfy the selected constraints.\n');
    return;
end

best_design = feasible_designs(1, :);

fprintf('\n=== HIGHEST-RANKED FEASIBLE DESIGN ===\n');
fprintf('Detector: %s (%d px width, %.1f um pixels)\n', ...
    best_design.det_name{1}, best_design.nx, best_design.pix_um);

fprintf('Spectrograph: R=%d, s1=%.0f um, rho=%.0f l/mm, alpha=%.1f deg\n', ...
    best_design.R, best_design.s1_um, ...
    best_design.rho_lpm, best_design.alpha_deg);

fprintf(' f1=%.1f mm, f2=%.1f mm, D1=%.1f mm, D2=%.1f mm\n', ...
    best_design.f1_mm, best_design.f2_mm, ...
    best_design.D1_mm, best_design.D2_mm);

fprintf(' Focal-length ratio=%.3f, slit-image magnification=%.3f\n', ...
    best_design.focal_length_ratio, best_design.slit_magnification);

detector_width_needed_mm = ...
    best_design.n_pixels_range * best_design.pix_um / 1000;

fprintf('\nPerformance: %.1f px/res element, %.0f pixels (%.1f mm width) needed\n', ...
    best_design.pixels_per_res, best_design.n_pixels_range, ...
    detector_width_needed_mm);

fprintf(' Uses %.0f/%.0f pixels (%.1f%% of detector)\n', ...
    best_design.n_pixels_range, best_design.nx, ...
    best_design.n_pixels_range / best_design.nx * 100);

linear_dispersion = best_design.dx_mm_per_nm;
fprintf(' Linear dispersion: %.3f mm/nm (%.2f nm/mm)\n', ...
    linear_dispersion, 1 / linear_dispersion);

%% Plot resolving power
if params.make_plot

    best_idx = best_design.config_id;

    figure();
    plot(results(best_idx).lambda_vec * 1e9, ...
        results(best_idx).R_all, 'b-', 'LineWidth', 2);
    hold on;

    plot([lambda_min, lambda_max] * 1e9, ...
        [best_design.R, best_design.R], ...
        'r--', 'LineWidth', 1.5);

    xlabel('wavelength [nm]');
    ylabel('resolving power');

    title(sprintf(['Resolving Power vs. Wavelength\n' ...
        'Detector: %s, R=%d, s1=%.1f um, rho=%.0f l/mm'], ...
        best_design.det_name{1}, best_design.R, ...
        best_design.s1_um, best_design.rho_lpm));

    grid on;
    box on;

    min_R = min(results(best_idx).R_all);
    max_R = max(results(best_idx).R_all);

    xlim([lambda_min, lambda_max] * 1e9);
    ylim([0.9 * min_R, 1.1 * max_R]);

    fprintf('\nResolving power analysis:\n');
    fprintf('  Minimum R across band: %.0f\n', min_R);
    fprintf('  Maximum R across band: %.0f\n', max_R);
end

end
