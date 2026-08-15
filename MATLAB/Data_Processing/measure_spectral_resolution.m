function resolution = measure_spectral_resolution(calibrated_file, varargin)
% MEASURE_SPECTRAL_RESOLUTION - Measure resolving power from an emission line
%
% Selects or accepts wavelength limits around an isolated calibrated line,
% estimates a local linear baseline, determines the FWHM by interpolation,
% and calculates R = lambda0/FWHM.
%
% Inputs:
%   calibrated_file  - CSV containing wavelength and intensity
%   'line_limits'    - [min,max] wavelength [nm]; interactive if empty
%   'make_plot'      - Display selection and result figures (default: true)
%
% Output:
%   resolution - Structure with centre wavelength, FWHM, R and sampling

p = inputParser;
addRequired(p, 'calibrated_file', @(x) ischar(x) || isstring(x));
addParameter(p, 'line_limits', [], ...
    @(x) isempty(x) || (isnumeric(x) && numel(x)==2 && x(2)>x(1)));
addParameter(p, 'make_plot', true, ...
    @(x) islogical(x) || ismember(x, [0,1]));

parse(p, calibrated_file, varargin{:});
params = p.Results;

calibrated_file = char(calibrated_file);
if ~isfile(calibrated_file)
    error('Calibrated spectrum not found: %s', calibrated_file);
end

data = readmatrix(calibrated_file);
if size(data,2) < 2
    error('Calibrated CSV must contain wavelength and intensity columns.');
end

% Accept either [wavelength,intensity] or [pixel,wavelength,intensity].
if size(data,2) >= 3
    wavelength = data(:,2);
    intensity = data(:,3);
else
    wavelength = data(:,1);
    intensity = data(:,2);
end

if isempty(params.line_limits)
    figure('Name', 'Resolution line selection');
    plot(wavelength, intensity, 'k', 'LineWidth', 1);
    xlabel('wavelength (nm)');
    ylabel('intensity');
    title('Select LEFT and RIGHT limits of one isolated emission line');
    grid on;

    disp('Click LEFT boundary');
    [x1, ~] = ginput(1);
    disp('Click RIGHT boundary');
    [x2, ~] = ginput(1);

    xmin = min(x1,x2);
    xmax = max(x1,x2);
else
    xmin = params.line_limits(1);
    xmax = params.line_limits(2);
end

idx = wavelength >= xmin & wavelength <= xmax;
wl = wavelength(idx);
signal = intensity(idx);

if numel(wl) < 7
    error('Selected line contains too few samples.');
end

n_edge = max(2, round(0.15*numel(wl)));
n_edge = min(n_edge, floor(numel(wl)/2));

wl_base = [wl(1:n_edge); wl(end-n_edge+1:end)];
signal_base = [signal(1:n_edge); signal(end-n_edge+1:end)];

baseline_coeff = polyfit(wl_base, signal_base, 1);
baseline = polyval(baseline_coeff, wl);
corrected = signal - baseline;

[peak_height, idx_peak] = max(corrected);
lambda0 = wl(idx_peak);
half_max = peak_height/2;

left_idx = find(corrected(1:idx_peak) <= half_max, 1, 'last');
right_rel = find(corrected(idx_peak:end) <= half_max, 1, 'first');

if isempty(left_idx) || isempty(right_rel)
    error('Could not determine both half-maximum crossings.');
end

right_idx = idx_peak + right_rel - 1;

if left_idx >= idx_peak || right_idx <= idx_peak
    error('Invalid half-maximum crossing geometry.');
end

lambda_left = interp1(...
    corrected(left_idx:left_idx+1), ...
    wl(left_idx:left_idx+1), ...
    half_max, 'linear');

lambda_right = interp1(...
    corrected(right_idx-1:right_idx), ...
    wl(right_idx-1:right_idx), ...
    half_max, 'linear');

fwhm_nm = lambda_right - lambda_left;
R_measured = lambda0 / fwhm_nm;

pixel_spacing_nm = median(abs(diff(wavelength)));
sampling_px = fwhm_nm / pixel_spacing_nm;

resolution = struct(...
    'lambda0_nm', lambda0, ...
    'fwhm_nm', fwhm_nm, ...
    'resolving_power', R_measured, ...
    'sampling_px_per_fwhm', sampling_px, ...
    'line_limits_nm', [xmin xmax], ...
    'baseline_coeff', baseline_coeff);

fprintf('Line centre: %.4f nm\n', lambda0);
fprintf('FWHM: %.4f nm\n', fwhm_nm);
fprintf('Measured resolving power: %.0f\n', R_measured);
fprintf('Sampling: %.2f px/FWHM\n', sampling_px);

if params.make_plot
    figure('Name', 'Measured spectral resolution');
    hold on;

    plot(wl, signal, 'k', 'LineWidth', 1.3);
    plot(wl, baseline, 'r--', 'LineWidth', 1.2);

    y_half_left = polyval(baseline_coeff, lambda_left) + half_max;
    y_half_right = polyval(baseline_coeff, lambda_right) + half_max;

    plot([lambda_left lambda_right], ...
         [y_half_left y_half_right], ...
         'm-', 'LineWidth', 2);

    xline(lambda0, 'k:');

    xlabel('wavelength (nm)');
    ylabel('intensity');
    title(sprintf('R = %.0f, FWHM = %.4f nm at %.2f nm', ...
        R_measured, fwhm_nm, lambda0));
    legend('Spectrum', 'Baseline', 'FWHM');
    grid on;
end

end
