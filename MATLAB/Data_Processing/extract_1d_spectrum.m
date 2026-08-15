function [spectrum_data, extraction_data] = extract_1d_spectrum(image_file, varargin)
% EXTRACT_1D_SPECTRUM - Extract a 1D spectrum from a 2D detector image
%
% Detects the illuminated spatial region, estimates a column-dependent
% background from a nearby dark strip, subtracts the background, and sums
% the corrected signal along the spatial direction.
%
% Inputs:
%   image_file             - Detector image path
%   'acquisition_mode'     - 'slit' or 'pinhole' (default: 'slit')
%   'subtract_dark'        - Apply dark-strip subtraction (default: true)
%   'slit_threshold'       - Relative slit-detection threshold (default: 0.20)
%   'slit_margin'          - Extra rows around detected slit (default: 5)
%   'pinhole_half_width'   - Half-width of pinhole extraction [rows] (default: 6)
%   'slit_dark_offset'     - Dark-strip offset for slit mode [rows] (default: 200)
%   'pinhole_dark_offset'  - Dark-strip offset for pinhole mode [rows] (default: 30)
%   'output_dir'           - Output directory (default: '1d_spectra')
%   'save_results'         - Save CSV, MAT and PNG outputs (default: true)
%   'make_plots'           - Display diagnostic figures (default: true)
%
% Outputs:
%   spectrum_data   - Table with detector pixel and extracted intensity
%   extraction_data - Structure with extraction regions and corrected data

p = inputParser;
addRequired(p, 'image_file', @(x) ischar(x) || isstring(x));
addParameter(p, 'acquisition_mode', 'slit', ...
    @(x) any(strcmpi(x, {'slit','pinhole'})));
addParameter(p, 'subtract_dark', true, ...
    @(x) islogical(x) || ismember(x, [0,1]));
addParameter(p, 'slit_threshold', 0.20, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
addParameter(p, 'slit_margin', 5, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'pinhole_half_width', 6, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'slit_dark_offset', 200, ...
    @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'pinhole_dark_offset', 30, ...
    @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'output_dir', '1d_spectra', ...
    @(x) ischar(x) || isstring(x));
addParameter(p, 'save_results', true, ...
    @(x) islogical(x) || ismember(x, [0,1]));
addParameter(p, 'make_plots', true, ...
    @(x) islogical(x) || ismember(x, [0,1]));

parse(p, image_file, varargin{:});
params = p.Results;

image_file = char(image_file);
if ~isfile(image_file)
    error('Detector image not found: %s', image_file);
end

raw = imread(image_file);
raw_d = double(raw);

if ndims(raw_d) ~= 2
    error('Expected a monochrome 2D detector image.');
end

[n_rows, n_cols] = size(raw_d);
row_profile = mean(raw_d, 2);

fprintf('Image size: %d rows x %d columns\n', n_rows, n_cols);

switch lower(params.acquisition_mode)

    case 'slit'
        threshold = median(row_profile) + params.slit_threshold * ...
            (max(row_profile) - median(row_profile));

        rows_above = find(row_profile > threshold);
        if isempty(rows_above)
            error('No slit signal detected above the selected threshold.');
        end

        row_top = max(1, min(rows_above) - params.slit_margin);
        row_bottom = min(n_rows, max(rows_above) + params.slit_margin);

        signal_rows = row_top:row_bottom;
        signal_height = numel(signal_rows);

        [dark_rows, dark_top, dark_bottom] = choose_dark_region(...
            row_top, row_bottom, signal_height, ...
            params.slit_dark_offset, n_rows);

        fprintf('Slit detected from row %d to %d (height = %d pixels)\n', ...
            row_top, row_bottom, signal_height);

    case 'pinhole'
        [~, row_c] = max(row_profile);

        row_top = max(1, row_c - params.pinhole_half_width);
        row_bottom = min(n_rows, row_c + params.pinhole_half_width);

        signal_rows = row_top:row_bottom;
        signal_height = numel(signal_rows);

        [dark_rows, dark_top, dark_bottom] = choose_dark_region(...
            row_top, row_bottom, signal_height, ...
            params.pinhole_dark_offset, n_rows);

        fprintf('Pinhole center at row %d, extraction rows %d-%d\n', ...
            row_c, row_top, row_bottom);
end

cropped = raw_d(signal_rows, :);
cropped_dark = raw_d(dark_rows, :);

if params.subtract_dark
    dark_spectrum = mean(cropped_dark, 1);
    cropped_corrected = cropped - dark_spectrum;
else
    dark_spectrum = zeros(1, n_cols);
    cropped_corrected = cropped;
end

spectrum_1d = sum(cropped_corrected, 1);
pixel = (1:n_cols)';

spectrum_data = table(pixel, spectrum_1d', ...
    'VariableNames', {'pixel','intensity'});

extraction_data = struct(...
    'raw_image_size', [n_rows, n_cols], ...
    'signal_rows', signal_rows, ...
    'dark_rows', dark_rows, ...
    'cropped', cropped, ...
    'cropped_dark', cropped_dark, ...
    'dark_spectrum', dark_spectrum, ...
    'cropped_corrected', cropped_corrected, ...
    'acquisition_mode', params.acquisition_mode, ...
    'subtract_dark', logical(params.subtract_dark));

fig_1d = [];
if params.make_plots
    figure('Name', 'Raw detector image');
    imshow(raw, []);
    title('Raw detector image');
    hold on;
    rectangle('Position', ...
        [1, row_top, n_cols, row_bottom-row_top+1], ...
        'EdgeColor', 'r', 'LineWidth', 1.5);
    rectangle('Position', ...
        [1, dark_top, n_cols, dark_bottom-dark_top+1], ...
        'EdgeColor', 'g', 'LineWidth', 1.5);
    hold off;

    figure('Name', 'Cropped signal strip');
    imshow(cropped, []);
    title('Signal region');

    figure('Name', '2D Spectral Map (corrected)');
    imagesc(cropped_corrected);
    xlabel('Spectral pixels');
    ylabel('Spatial pixels (rows)');
    title('2D Spectral Map (background subtracted)');
    colorbar;

    n_rows_cropped = size(cropped_corrected, 1);
    if n_rows_cropped <= 30
        step = 1;
    else
        step = max(1, round(n_rows_cropped/10));
    end

    tick_idx = 1:step:n_rows_cropped;
    yticks(tick_idx);
    yticklabels(signal_rows(tick_idx));

    fig_1d = figure('Name', '1D Spectrum');
    plot(pixel, spectrum_1d, 'b-', 'LineWidth', 1.5);
    xlim([1 n_cols]);
    xlabel('Spectral pixels');
    ylabel('Intensity (counts)');
    title('1D Spectrum');
    grid on;
end

if params.save_results
    output_dir = char(params.output_dir);
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    [~, image_base, ~] = fileparts(image_file);

    writetable(spectrum_data, ...
        fullfile(output_dir, [image_base '_1d_spectrum.csv']));

    save(fullfile(output_dir, [image_base '_extraction.mat']), ...
        'spectrum_data', 'extraction_data');

    if ~isempty(fig_1d)
        exportgraphics(fig_1d, ...
            fullfile(output_dir, [image_base '_1d_spectrum.png']), ...
            'Resolution', 300);
    end
end

end

function [dark_rows, dark_top, dark_bottom] = choose_dark_region(...
    row_top, row_bottom, signal_height, offset, n_rows)

dark_top = row_top + offset;
dark_bottom = dark_top + signal_height - 1;

if dark_bottom > n_rows
    dark_bottom = row_top - offset;
    dark_top = dark_bottom - signal_height + 1;
end

if dark_top < 1 || dark_bottom > n_rows
    error('Could not place a dark region of equal height inside the detector.');
end

dark_rows = dark_top:dark_bottom;

end
