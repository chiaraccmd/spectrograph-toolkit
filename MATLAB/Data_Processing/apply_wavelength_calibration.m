function calibrated_spectrum = apply_wavelength_calibration(spectrum_file, calibration_file, varargin)
% APPLY_WAVELENGTH_CALIBRATION - Apply a saved wavelength solution
%
% Applies polynomial calibration coefficients to a two-column extracted
% spectrum and optionally saves the calibrated spectrum and figure.
%
% Inputs:
%   spectrum_file    - Two-column CSV: pixel, intensity
%   calibration_file - MAT file produced by wavelength_calibration
%   'output_dir'     - Output directory (default: 'calibrated_spectra')
%   'output_name'    - Base output name; defaults to input base name
%   'save_results'   - Save CSV and PNG (default: true)
%   'make_plot'      - Display calibrated spectrum (default: true)
%
% Output:
%   calibrated_spectrum - Table with pixel, wavelength_nm and intensity

p = inputParser;
addRequired(p, 'spectrum_file', @(x) ischar(x) || isstring(x));
addRequired(p, 'calibration_file', @(x) ischar(x) || isstring(x));
addParameter(p, 'output_dir', 'calibrated_spectra', ...
    @(x) ischar(x) || isstring(x));
addParameter(p, 'output_name', '', ...
    @(x) ischar(x) || isstring(x));
addParameter(p, 'save_results', true, ...
    @(x) islogical(x) || ismember(x, [0,1]));
addParameter(p, 'make_plot', true, ...
    @(x) islogical(x) || ismember(x, [0,1]));

parse(p, spectrum_file, calibration_file, varargin{:});
params = p.Results;

spectrum_file = char(spectrum_file);
calibration_file = char(calibration_file);

if ~isfile(spectrum_file)
    error('Spectrum file not found: %s', spectrum_file);
end

if ~isfile(calibration_file)
    error('Calibration file not found: %s', calibration_file);
end

data = readmatrix(spectrum_file);
if size(data,2) < 2
    error('Spectrum CSV must contain at least two columns: pixel and intensity.');
end

cal = load(calibration_file);

if isfield(cal, 'calibration')
    calibration = cal.calibration;
    if ~isfield(calibration, 'p')
        error('Calibration structure does not contain polynomial coefficients.');
    end
    p_fit = calibration.p;
elseif isfield(cal, 'p')
    p_fit = cal.p;
else
    error('Calibration file does not contain polynomial coefficients.');
end

pixel = data(:,1);
intensity = data(:,2);
wavelength_nm = polyval(p_fit, pixel);

calibrated_spectrum = table(pixel, wavelength_nm, intensity);

fig = [];
if params.make_plot
    fig = figure('Name', 'Calibrated spectrum');
    plot(wavelength_nm, intensity, 'k', 'LineWidth', 1);
    xlabel('wavelength (nm)');
    ylabel('intensity');
    title('Calibrated spectrum');
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

    writetable(calibrated_spectrum, ...
        fullfile(output_dir, [base_name '_calibrated.csv']));

    if ~isempty(fig)
        exportgraphics(fig, ...
            fullfile(output_dir, [base_name '_calibrated_spectrum.png']), ...
            'Resolution', 300);
    end
end

end
