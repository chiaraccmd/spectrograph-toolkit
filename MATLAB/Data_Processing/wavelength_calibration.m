function [calibration, calibrated_spectrum] = wavelength_calibration(spectrum_file, varargin)
% WAVELENGTH_CALIBRATION - Fit a detector-pixel wavelength calibration
%
% Detects emission peaks in a Neon reference spectrum, matches them to
% reference wavelengths, fits a polynomial wavelength solution, and reports
% calibration residuals in wavelength and detector pixels.
%
% Inputs:
%   spectrum_file         - Two-column CSV: pixel, intensity
%   'poly_order'          - Polynomial order (default: 2)
%   'wavelength_range'    - Calibration range [nm] (default: [560 721])
%   'peak_use_idx'        - Optional detected-peak indices for manual matching
%   'ref_use_idx'         - Optional reference-line indices for manual matching
%   'output_dir'          - Output directory (default: 'calibrated_spectra')
%   'output_name'         - Base output name; defaults to input base name
%   'save_results'        - Save calibration MAT/CSV/PNG (default: true)
%   'make_plots'          - Display calibration diagnostics (default: true)
%
% Outputs:
%   calibration          - Polynomial coefficients and residual statistics
%   calibrated_spectrum  - Table with pixel, wavelength and intensity

p = inputParser;
addRequired(p, 'spectrum_file', @(x) ischar(x) || isstring(x));
addParameter(p, 'poly_order', 2, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'wavelength_range', [560 721], ...
    @(x) isnumeric(x) && numel(x) == 2 && x(2) > x(1));
addParameter(p, 'peak_use_idx', [], @isnumeric);
addParameter(p, 'ref_use_idx', [], @isnumeric);
addParameter(p, 'output_dir', 'calibrated_spectra', ...
    @(x) ischar(x) || isstring(x));
addParameter(p, 'output_name', '', ...
    @(x) ischar(x) || isstring(x));
addParameter(p, 'save_results', true, ...
    @(x) islogical(x) || ismember(x, [0,1]));
addParameter(p, 'make_plots', true, ...
    @(x) islogical(x) || ismember(x, [0,1]));

parse(p, spectrum_file, varargin{:});
params = p.Results;

spectrum_file = char(spectrum_file);
if ~isfile(spectrum_file)
    error('Spectrum file not found: %s', spectrum_file);
end

data = readmatrix(spectrum_file);
if size(data,2) < 2
    error('Spectrum CSV must contain at least two columns: pixel and intensity.');
end

pixels = data(:,1)';
neon_spec = data(:,2)';

ne_ref_all = [...
    585.25, 588.19, 594.48, 597.55, 603.00, 607.43, ...
    609.62, 614.31, 616.36, 621.73, 626.65, 630.48, ...
    633.44, 638.30, 640.23, 650.65, 653.29, 659.90, ...
    667.83, 671.70, 692.95];

wl_min = params.wavelength_range(1);
wl_max = params.wavelength_range(2);
ne_ref = ne_ref_all(ne_ref_all >= wl_min & ne_ref_all <= wl_max);

threshold = median(neon_spec) + 0.3 * std(neon_spec);

[~, peak_pix] = findpeaks(neon_spec, pixels, ...
    'MinPeakHeight', threshold, ...
    'MinPeakDistance', 3);

fprintf('Found %d peaks and %d reference lines in [%.0f-%.0f] nm\n', ...
    numel(peak_pix), numel(ne_ref), wl_min, wl_max);

if ~isempty(params.peak_use_idx) || ~isempty(params.ref_use_idx)
    if numel(params.peak_use_idx) ~= numel(params.ref_use_idx)
        error('peak_use_idx and ref_use_idx must contain the same number of entries.');
    end

    matched_pix = peak_pix(params.peak_use_idx);
    matched_wl = ne_ref(params.ref_use_idx);

elseif numel(peak_pix) == numel(ne_ref)
    matched_pix = peak_pix;
    matched_wl = ne_ref;

else
    error(['Detected peaks (%d) and reference lines (%d) do not match. ' ...
           'Inspect the spectrum and supply peak_use_idx and ref_use_idx.'], ...
          numel(peak_pix), numel(ne_ref));
end

if numel(matched_pix) < params.poly_order + 1
    error('At least %d matched lines are required for polynomial order %d.', ...
        params.poly_order + 1, params.poly_order);
end

matched_intensity = interp1(pixels, neon_spec, matched_pix, 'linear');

p_fit = polyfit(matched_pix, matched_wl, params.poly_order);
wl_axis = polyval(p_fit, pixels);

res_nm = matched_wl - polyval(p_fit, matched_pix);
rms_res_nm = sqrt(mean(res_nm.^2));

dp = polyder(p_fit);
local_dispersion = polyval(dp, matched_pix);
res_pix = res_nm ./ local_dispersion;
rms_res_pix = sqrt(mean(res_pix.^2));

detector_centre = mean([min(pixels), max(pixels)]);
dispersion_centre = polyval(dp, detector_centre);

fprintf('Matched %d lines\n', numel(matched_pix));
fprintf('Calibration RMS: %.4f nm\n', rms_res_nm);
fprintf('Calibration RMS: %.4f pixel\n', rms_res_pix);
fprintf('Maximum residual: %.4f nm\n', max(abs(res_nm)));
fprintf('Local dispersion at detector centre: %.4f nm/px\n', dispersion_centre);

calibration = struct(...
    'p', p_fit, ...
    'poly_order', params.poly_order, ...
    'matched_pix', matched_pix, ...
    'matched_wl', matched_wl, ...
    'residual_nm', res_nm, ...
    'residual_pix', res_pix, ...
    'rms_res_nm', rms_res_nm, ...
    'rms_res_pix', rms_res_pix, ...
    'dispersion_centre_nm_per_px', dispersion_centre);

calibrated_spectrum = table(...
    pixels', wl_axis', neon_spec', ...
    'VariableNames', {'pixel','wavelength_nm','intensity'});

fig_cal = [];
if params.make_plots
    figure('Name', 'Peak identification');
    plot(pixels, neon_spec, 'k', 'LineWidth', 1);
    hold on;
    plot(matched_pix, matched_intensity, 'rv', ...
        'MarkerSize', 5, 'MarkerFaceColor', 'r');

    for k = 1:numel(matched_pix)
        text(matched_pix(k), matched_intensity(k)*1.06, ...
            sprintf('%.2f', matched_wl(k)), ...
            'FontSize', 7, ...
            'HorizontalAlignment', 'center', ...
            'Color', 'r');
    end

    xlabel('pixels');
    ylabel('intensity');
    title('Matched Neon lines');
    grid on;

    figure('Name', 'Calibration curve');
    plot(matched_pix, matched_wl, 'bo', ...
        'MarkerSize', 8, 'MarkerFaceColor', 'b', ...
        'DisplayName', 'Matched lines');
    hold on;
    plot(pixels, wl_axis, 'r-', ...
        'LineWidth', 1.5, ...
        'DisplayName', sprintf('Polynomial fit, order %d', params.poly_order));
    xlabel('pixels');
    ylabel('wavelength (nm)');
    title(sprintf('Calibration (RMS = %.4f nm)', rms_res_nm));
    legend('Location', 'best');
    grid on;

    figure('Name', 'Residuals');
    scatter(matched_pix, res_nm, 20, 'filled');
    hold on;
    yline(0, 'r--', 'LineWidth', 1.5);
    xlabel('pixel');
    ylabel('residual (nm)');
    title('Calibration residuals');
    grid on;

    fig_cal = figure('Name', 'Calibrated spectrum');
    plot(wl_axis, neon_spec, 'k', 'LineWidth', 1);
    xlabel('wavelength (nm)');
    ylabel('intensity');
    title('Calibrated Neon spectrum');
    grid on;
end

if params.save_results
    output_dir = char(params.output_dir);
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    [~, input_base, ~] = fileparts(spectrum_file);
    if strlength(string(params.output_name)) == 0
        base_name = erase(input_base, '_1d_spectrum');
    else
        base_name = char(params.output_name);
    end

    writetable(calibrated_spectrum(:, {'pixel','wavelength_nm'}), ...
        fullfile(output_dir, [base_name '_wavelength_axis.csv']));

    save(fullfile(output_dir, [base_name '_calibration.mat']), ...
        'calibration');

    if ~isempty(fig_cal)
        exportgraphics(fig_cal, ...
            fullfile(output_dir, [base_name '_calibrated_spectrum.png']), ...
            'Resolution', 300);
    end
end

end
