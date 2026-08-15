function peak = integrate_spectral_peak(spectrum_file, varargin)
% INTEGRATE_SPECTRAL_PEAK - Integrate a baseline-corrected spectral peak
%
% Integrates a selected spectral feature after fitting a local linear
% baseline to samples near the two edges of the selected interval.
%
% Inputs:
%   spectrum_file    - CSV containing x-axis and intensity
%   'limits'         - [min,max] integration limits; interactive if empty
%   'x_type'         - 'auto', 'pixel', or 'wavelength' (default: 'auto')
%   'grating_lines'  - Optional grating density [lines/mm] for plot title
%   'make_plot'      - Display selection and integration figures
%
% Output:
%   peak - Structure with position, height, area, limits and baseline

p = inputParser;
addRequired(p, 'spectrum_file', @(x) ischar(x) || isstring(x));
addParameter(p, 'limits', [], ...
    @(x) isempty(x) || (isnumeric(x) && numel(x)==2 && x(2)>x(1)));
addParameter(p, 'x_type', 'auto', ...
    @(x) any(strcmpi(x, {'auto','pixel','wavelength'})));
addParameter(p, 'grating_lines', [], ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'make_plot', true, ...
    @(x) islogical(x) || ismember(x, [0,1]));

parse(p, spectrum_file, varargin{:});
params = p.Results;

spectrum_file = char(spectrum_file);
if ~isfile(spectrum_file)
    error('Spectrum file not found: %s', spectrum_file);
end

data = readmatrix(spectrum_file);
if size(data,2) < 2
    error('Spectrum CSV must contain at least two columns.');
end

switch lower(params.x_type)
    case 'wavelength'
        if size(data,2) >= 3
            x = data(:,2)';
            counts = data(:,3)';
        else
            x = data(:,1)';
            counts = data(:,2)';
        end
        x_label = 'wavelength (nm)';
        area_units = 'counts*nm';

    case 'pixel'
        x = data(:,1)';
        counts = data(:,end)';
        x_label = 'pixels';
        area_units = 'counts*pixel';

    case 'auto'
        if size(data,2) >= 3
            x = data(:,2)';
            counts = data(:,3)';
            x_label = 'wavelength (nm)';
            area_units = 'counts*nm';
        else
            x = data(:,1)';
            counts = data(:,2)';
            x_label = 'x';
            area_units = 'counts*x-unit';
        end
end

if isempty(params.limits)
    figure('Name', 'Peak selection');
    plot(x, counts, 'b', 'LineWidth', 1.5);
    xlabel(x_label);
    ylabel('counts');
    title('Select LEFT and RIGHT limits of the peak');
    grid on;

    disp('Click LEFT boundary');
    [x1, ~] = ginput(1);
    disp('Click RIGHT boundary');
    [x2, ~] = ginput(1);

    xmin = min(x1,x2);
    xmax = max(x1,x2);
else
    xmin = params.limits(1);
    xmax = params.limits(2);
end

idx = x >= xmin & x <= xmax;
xp = x(idx);
yp = counts(idx);

if numel(xp) < 5
    error('Selected peak region contains too few samples.');
end

n_edge = max(2, round(0.1*numel(xp)));
n_edge = min(n_edge, floor(numel(xp)/2));

x_base = [xp(1:n_edge), xp(end-n_edge+1:end)];
y_base = [yp(1:n_edge), yp(end-n_edge+1:end)];

baseline_coeff = polyfit(x_base, y_base, 1);
baseline = polyval(baseline_coeff, xp);
ycorr = yp - baseline;

[peak_height, idx_max] = max(ycorr);
peak_position = xp(idx_max);
peak_area = trapz(xp, ycorr);

peak = struct(...
    'position', peak_position, ...
    'height_counts', peak_height, ...
    'area', peak_area, ...
    'area_units', area_units, ...
    'limits', [xmin xmax], ...
    'baseline_coeff', baseline_coeff);

fprintf('Peak position = %.4f\n', peak_position);
fprintf('Peak height = %.3f counts\n', peak_height);
fprintf('Peak integration limits: %.4f - %.4f\n', xmin, xmax);
fprintf('Integrated area = %.3f %s\n', peak_area, area_units);

if params.make_plot
    figure('Name', 'Integrated peak');
    hold on;

    plot(xp, yp, 'b', 'LineWidth', 1.5);
    plot(xp, baseline, 'r--', 'LineWidth', 2);

    fill([xp fliplr(xp)], ...
         [yp fliplr(baseline)], ...
         'c', 'EdgeColor', 'none', 'FaceAlpha', 0.3);

    xlabel(x_label);
    ylabel('counts');

    if isempty(params.grating_lines)
        title(sprintf('Integrated Area = %.2f', peak_area));
    else
        title(sprintf('Integrated Area = %.2f (%.0f l/mm)', ...
            peak_area, params.grating_lines));
    end

    legend('Spectrum', 'Baseline', 'Integrated peak');
    grid on;
end

end
