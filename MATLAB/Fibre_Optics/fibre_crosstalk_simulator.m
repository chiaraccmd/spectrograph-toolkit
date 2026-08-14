function [crosstalk_results, analysis_data] = fibre_crosstalk_simulator(model_type, varargin)
% FIBRE_CROSSTALK_SIMULATOR - Comprehensive fibre crosstalk analysis toolkit
%
% Analyzes crosstalk between optical fibres using multiple PSF models:
% - Airy disk model (diffraction-limited)
% - Gaussian beam model (modal propagation)
% - Dispersed spectrum model (multi-wavelength)
% - Pixel-integrated analysis (detector effects)
%
% Physical principles:
% - Fourier optics and diffraction theory
% - Gaussian beam propagation
% - Grating dispersion and spectral imaging
% - Detector sampling and pixel integration
%
% Inputs:
%   model_type - 'airy', 'gaussian', 'dispersed', or 'all'
%   Optional parameter/value pairs override defaults
%
% Outputs:
%   crosstalk_results - Structure with crosstalk metrics for all models
%   analysis_data     - Detailed simulation data and PSF arrays
%
% Example:
%   % Airy model analysis
%   [results, data] = fibre_crosstalk_simulator('airy', 'fibre_separation', 25e-6);
%
%   % Complete analysis suite
%   [results, data] = fibre_crosstalk_simulator('all', 'wavelength', 1.55e-6);

%% Input parsing and parameter initialization
% =========================================================================
p = inputParser;
addRequired(p, 'model_type', @(x) any(validatestring(x, {'airy', 'gaussian', 'dispersed', 'all'})));
addParameter(p, 'fibre_separation', 25e-6, @isnumeric);           % Core separation [m]
addParameter(p, 'fibre_diameter', 7.3e-6, @isnumeric);            % Core diameter [m]
addParameter(p, 'wavelength', 1.55e-6, @isnumeric);               % Analysis wavelength [m]
addParameter(p, 'pixel_size', 18e-6, @isnumeric);                 % Detector pixel pitch [m]
addParameter(p, 'nPix', 2000, @isnumeric);                        % Detector pixel count
addParameter(p, 'F1', 4.55, @isnumeric);                          % Collimator f-number
addParameter(p, 'F2', 19.21, @isnumeric);                         % Camera f-number
addParameter(p, 'anamorphism', 1, @isnumeric);                    % Anamorphic factor
addParameter(p, 'grating_density', 412e3, @isnumeric);            % Grating lines/m
addParameter(p, 'incidence_angle', deg2rad(19.7056), @isnumeric); % Incidence angle [rad]
addParameter(p, 'camera_focal_length', 0.2722, @isnumeric);       % Camera focal length [m]

% Wavelength bands
addParameter(p, 'L_min_Y', 0.960e-6, @isnumeric);
addParameter(p, 'L_max_Y', 1.115e-6, @isnumeric);
addParameter(p, 'L_min_J', 1.150e-6, @isnumeric);
addParameter(p, 'L_max_J', 1.345e-6, @isnumeric);
addParameter(p, 'L_min_H', 1.490e-6, @isnumeric);
addParameter(p, 'L_max_H', 1.780e-6, @isnumeric);

% Gaussian model parameters
addParameter(p, 'MFD_1310', 9.2e-6, @isnumeric);                  % Mode field diameter [m]
addParameter(p, 'lambda_ref', 1310e-9, @isnumeric);               % Reference wavelength [m]

% Simulation parameters
addParameter(p, 'grid_size', 2048, @isnumeric);                   % Simulation grid size
addParameter(p, 'subsampling_factor', 10, @isnumeric);            % Pixel subsampling
addParameter(p, 'fov_pixels', 20, @isnumeric);                    % Field of view [pixels]

parse(p, model_type, varargin{:});
params = p.Results;

% Calculate derived parameters
params.M = params.F2 / params.F1;                                 % System magnification
params.fibre_radius_detector = (params.fibre_diameter/2) * params.M; % Image plane radius

fprintf('=== Fibre Crosstalk Analysis ===\n');
fprintf('Model: %s, Separation: %.1f μm, λ: %.1f nm\n', ...
    model_type, params.fibre_separation*1e6, params.wavelength*1e9);

%% Model Selection and Execution
% =========================================================================
crosstalk_results = struct();
analysis_data = struct();

switch lower(model_type)
    case 'airy'
        [crosstalk_results.airy, analysis_data.airy] = airy_crosstalk_analysis(params);
        
    case 'gaussian'
        [crosstalk_results.gaussian, analysis_data.gaussian] = gaussian_crosstalk_analysis(params);
        
    case 'dispersed'
        [crosstalk_results.dispersed, analysis_data.dispersed] = dispersed_crosstalk_analysis(params);
        
    case 'all'
        % Run all models for comprehensive comparison
        [crosstalk_results.airy, analysis_data.airy] = airy_crosstalk_analysis(params);
        [crosstalk_results.gaussian, analysis_data.gaussian] = gaussian_crosstalk_analysis(params);
        [crosstalk_results.dispersed, analysis_data.dispersed] = dispersed_crosstalk_analysis(params);
        
        % Comparative analysis
        crosstalk_results.comparison = compare_crosstalk_models(crosstalk_results, params);
end

fprintf('\n=== Analysis Complete ===\n');

end

%% Airy PSF Crosstalk Analysis
% =========================================================================
function [results, data] = airy_crosstalk_analysis(params)
% AIRY_CROSSTALK_ANALYSIS - Diffraction-limited crosstalk using Airy disk model

    fprintf('\n--- Airy PSF Crosstalk Analysis ---\n');
    
    %% Continuous PSF Analysis
    % =====================================================================
    % Airy PSF function definition
    airy_psf_continuous = @(u) airy_intensity(u);
    
    % Grid setup for continuous analysis
    fov = 100e-6 * params.M; % Field of view [m]
    dx = fov / params.grid_size;
    x = linspace(-fov/2, fov/2, params.grid_size);
    [X, Y] = meshgrid(x, x);
    
    % Apply anamorphism to radial coordinate
    Rgrid = sqrt(X.^2 + (Y/params.anamorphism).^2);
    
    % Airy first-zero radius at the requested analysis wavelength
    r0 = 1.22 * params.F2 * params.wavelength;
    fprintf('Airy radius: r₀ = %.2f μm\n', r0*1e6);
    
    % Compute and normalize PSF
    I_continuous = airy_psf_continuous(Rgrid / r0);
    I_continuous = I_continuous / sum(I_continuous(:));
    
    %% Crosstalk vs Distance Analysis
    % =====================================================================
    d_max = max(30e-6, 1.2 * params.fibre_separation);
    d_steps = linspace(0, d_max, 100); % Separation sweep [m]
    d_detector = d_steps * params.M; % Image plane separation
    
    crosstalk_continuous = zeros(size(d_detector));
    for k = 1:length(d_detector)
        shifted_mask = ((X - d_detector(k))/r0).^2 + ...
                       (Y/(r0*params.anamorphism)).^2 <= 1;
        crosstalk_continuous(k) = sum(I_continuous(shifted_mask));
    end
    
    %% Pixel-Integrated Analysis
    % =====================================================================
    fprintf('Performing pixel-integrated analysis...\n');
    
    % Create function handle for pixel integration
    airy_psf_pixel = @(u) airy_intensity(u);
    
    [crosstalk_pixel, pixel_d_steps, I_pixel, x_pixel, y_pixel] = pixel_integrated_analysis(...
        params, airy_psf_pixel, 'airy', r0);
    
    %% Key Distance Analysis
    % =====================================================================
    key_separation = params.fibre_separation * params.M;
    [~, idx] = min(abs(d_detector - key_separation));
    crosstalk_at_key = crosstalk_continuous(idx) * 100;
    
    fprintf('Crosstalk at %.1f μm separation: %.3f%% (continuous)\n', ...
        params.fibre_separation*1e6, crosstalk_at_key);
    
    %% Visualization
    % =====================================================================
    if params.fibre_separation <= d_max
        fig_airy = visualize_airy_analysis(params, I_continuous, x, ...
            crosstalk_continuous, d_steps, key_separation, crosstalk_at_key, r0);
        data.figures.airy_continuous = fig_airy;
    end
    
    %% Prepare Outputs
    % =====================================================================
    results.continuous = struct(...
        'separations', d_steps, ...
        'crosstalk_percent', crosstalk_continuous * 100, ...
        'airy_radius', r0, ...
        'key_separation_crosstalk', crosstalk_at_key);
    
    results.pixel_integrated = struct(...
        'separations', pixel_d_steps, ...
        'crosstalk_percent', crosstalk_pixel * 100);
    
    data.psf_continuous = I_continuous;
    data.psf_pixel = I_pixel;
    data.grid = struct('x_continuous', x, 'x_pixel', x_pixel, 'y_pixel', y_pixel);
    data.parameters = params;
end

%% Gaussian PSF Crosstalk Analysis
% =========================================================================
function [results, data] = gaussian_crosstalk_analysis(params)
% GAUSSIAN_CROSSTALK_ANALYSIS - Gaussian beam propagation crosstalk model

    fprintf('\n--- Gaussian PSF Crosstalk Analysis ---\n');
    
    % Gaussian PSF function
    gaussian_psf = @(X, Y, wx, wy) exp(-2*(X.^2/wx^2 + Y.^2/wy^2));
    
    % Wavelength-dependent beam waist calculation
    calc_waist = @(lambda) (params.MFD_1310/2) * params.M * (lambda/params.lambda_ref);
    
    %% Grid Setup
    fov = 100e-6 * params.M;
    x = linspace(-fov/2, fov/2, params.grid_size);
    [X, Y] = meshgrid(x, x);
    
    %% Beam Parameter Calculation
    wx = calc_waist(params.wavelength);                      % Spatial direction
    wy = calc_waist(params.wavelength) * params.anamorphism; % Spectral direction
    
    fprintf('Gaussian waists: wx = %.2f μm, wy = %.2f μm\n', wx*1e6, wy*1e6);
    
    %% PSF Generation and Normalization
    I_gaussian = gaussian_psf(X, Y, wx, wy);
    I_gaussian = I_gaussian / sum(I_gaussian(:));
    
    %% Crosstalk vs Distance
    d_max = max(30e-6, 1.2 * params.fibre_separation);
    d_steps = linspace(0, d_max, 50);
    d_detector = d_steps * params.M;
    
    crosstalk_gaussian = zeros(size(d_detector));
    for k = 1:length(d_detector)
        shifted_mask = ((X - d_detector(k)).^2/(wx^2) + ...
                       Y.^2/(wy^2)) <= 1;
        crosstalk_gaussian(k) = sum(I_gaussian(shifted_mask));
    end
    
    %% Key Separation Analysis
    key_separation = params.fibre_separation * params.M;
    [~, idx] = min(abs(d_detector - key_separation));
    crosstalk_at_key = crosstalk_gaussian(idx) * 100;
    
    fprintf('Crosstalk at %.1f μm separation: %.3f%% (Gaussian)\n', ...
        params.fibre_separation*1e6, crosstalk_at_key);
    
    %% Visualization
    if params.fibre_separation <= d_max
        fig_gaussian = visualize_gaussian_analysis(params, I_gaussian, x, ...
            crosstalk_gaussian, d_steps, key_separation, crosstalk_at_key, wx, wy);
        data.figures.gaussian = fig_gaussian;
    end
    
    %% Prepare Outputs
    results = struct(...
        'separations', d_steps, ...
        'crosstalk_percent', crosstalk_gaussian * 100, ...
        'beam_waists', struct('wx', wx, 'wy', wy), ...
        'key_separation_crosstalk', crosstalk_at_key);
    
    data.psf_gaussian = I_gaussian;
    data.grid = x;
    data.parameters = params;
end

%% Dispersed Spectrum Crosstalk Analysis
% =========================================================================
function [results, data] = dispersed_crosstalk_analysis(params)
% DISPERSED_CROSSTALK_ANALYSIS - Multi-wavelength crosstalk with dispersion

    fprintf('\n--- Dispersed Spectrum Crosstalk Analysis ---\n');
    
    % Airy PSF function
    airy_psf = @(u) airy_intensity(u);
    
    %% Spectral Dispersion Calculation
    % =====================================================================
    lambda_samples = linspace(params.L_min_H, params.L_max_H, 15);
    
    % Diffraction angles and spectral positions
    beta = asin(params.grating_density * lambda_samples - sin(params.incidence_angle));
    beta_0 = asin(params.grating_density * mean(lambda_samples) - sin(params.incidence_angle));
    y_pos = params.camera_focal_length * (tan(beta) - tan(beta_0)); % Spectral positions
    
    fprintf('Spectral range: %.1f-%.1f nm (%d samples)\n', ...
        min(lambda_samples)*1e9, max(lambda_samples)*1e9, length(lambda_samples));
    
    %% Multi-Wavelength PSF Stack
    % =====================================================================
    fov = 100e-6 * params.M;
    x = linspace(-fov/2, fov/2, params.grid_size);
    y = linspace(min(y_pos)-10e-6, max(y_pos)+10e-6, params.grid_size);
    [X, Y] = meshgrid(x, y);
    
    PSF_stack = zeros(params.grid_size, params.grid_size, length(lambda_samples));
    r0_values = zeros(size(lambda_samples));
    
    for i = 1:length(lambda_samples)
        r0_values(i) = 1.22 * params.F2 * lambda_samples(i);
        R = sqrt(X.^2 + ((Y - y_pos(i))/params.anamorphism).^2);
        PSF_stack(:,:,i) = airy_psf(R / r0_values(i)); % Normalized argument
        PSF_stack(:,:,i) = PSF_stack(:,:,i) / sum(PSF_stack(:,:,i), 'all');
    end
    
    %% Combined PSF and Crosstalk Analysis
    % =====================================================================
    combined_psf = sum(PSF_stack, 3);
    combined_psf = combined_psf / sum(combined_psf(:));
    
    d_max = max(30e-6, 1.2 * params.fibre_separation);
    d_steps = linspace(0, d_max, 20);
    d_detector = d_steps * params.M;
    crosstalk_dispersed = zeros(size(d_detector));
    
    for k = 1:length(d_detector)
        total_flux = 0;
        leaked_flux = 0;
        
        for i = 1:length(lambda_samples)
            R_original = sqrt(X.^2 + ((Y - y_pos(i))/params.anamorphism).^2);
            
            R_shifted = sqrt((X - d_detector(k)).^2 + ...
                             ((Y - y_pos(i))/params.anamorphism).^2);
            mask_shifted = R_shifted <= r0_values(i);
            
            original_flux = sum(PSF_stack(:,:,i), 'all');
            leaked_flux = leaked_flux + sum(PSF_stack(:,:,i).*mask_shifted, 'all');
            total_flux = total_flux + original_flux;
        end
        
        crosstalk_dispersed(k) = leaked_flux / total_flux;
    end
    
    %% Pixel-Integrated Analysis
    % =====================================================================
    [crosstalk_pixel, pixel_d_steps, I_pixel, x_pixel, y_pixel] = pixel_integrated_dispersed_analysis(...
        params, lambda_samples, y_pos, r0_values);
    
    %% Key Separation Results
    key_separation = params.fibre_separation * params.M;
    [~, idx] = min(abs(d_detector - key_separation));
    crosstalk_at_key = crosstalk_dispersed(idx) * 100;
    
    fprintf('Crosstalk at %.1f μm separation: %.3f%% (dispersed)\n', ...
        params.fibre_separation*1e6, crosstalk_at_key);
    
    %% Visualization
    fig_dispersed = visualize_dispersed_analysis(params, combined_psf, x, y, ...
        y_pos, r0_values, lambda_samples, key_separation, crosstalk_at_key, ...
        d_steps, crosstalk_dispersed);
    data.figures.dispersed = fig_dispersed;
    
    %% Prepare Outputs
    results.continuous = struct(...
        'separations', d_steps, ...
        'crosstalk_percent', crosstalk_dispersed * 100, ...
        'key_separation_crosstalk', crosstalk_at_key);
    
    results.pixel_integrated = struct(...
        'separations', pixel_d_steps, ...
        'crosstalk_percent', crosstalk_pixel * 100);
    
    data.psf_stack = PSF_stack;
    data.combined_psf = combined_psf;
    data.spectral_data = struct(...
        'wavelengths', lambda_samples, ...
        'spectral_positions', y_pos, ...
        'airy_radii', r0_values);
    data.parameters = params;
end

%% Pixel Integration Core Function
% =========================================================================
function [crosstalk_pixel, d_steps, I_pixel, x_pixel, y_pixel] = pixel_integrated_analysis(params, psf_function, model_type, r0)
% PIXEL_INTEGRATED_ANALYSIS - Account for detector pixelation effects

    fov = params.fov_pixels * params.pixel_size;
    n_sub = params.fov_pixels * params.subsampling_factor;
    d_sub = params.pixel_size / params.subsampling_factor;

    x_subsampled = ((0:n_sub-1) + 0.5) * d_sub - fov/2;
    [X_subs, Y_subs] = meshgrid(x_subsampled, x_subsampled);

    switch model_type
        case 'airy'
            R_subs = sqrt(X_subs.^2 + (Y_subs/params.anamorphism).^2) / r0;
            I_subsampled = psf_function(R_subs);
        case 'gaussian'
            calc_waist = @(lambda) (params.MFD_1310/2) * params.M * (lambda/params.lambda_ref);
            wx = calc_waist(params.wavelength);
            wy = calc_waist(params.wavelength) * params.anamorphism;
            I_subsampled = exp(-2*(X_subs.^2/wx^2 + Y_subs.^2/wy^2));
    end

    I_subsampled = I_subsampled / sum(I_subsampled(:));

    I_pixel = zeros(params.fov_pixels, params.fov_pixels);
    q = params.subsampling_factor;
    for i = 1:params.fov_pixels
        for j = 1:params.fov_pixels
            idx_x = (i-1)*q + 1 : i*q;
            idx_y = (j-1)*q + 1 : j*q;
            I_pixel(j,i) = sum(I_subsampled(idx_y, idx_x), 'all');
        end
    end
    I_pixel = I_pixel / sum(I_pixel(:));

    x_pixel = ((0:params.fov_pixels-1) + 0.5) * params.pixel_size - fov/2;
    y_pixel = x_pixel;

    d_max = max(30e-6, 1.2 * params.fibre_separation);
    d_steps = linspace(0, d_max, 30);
    d_detector = d_steps * params.M;
    crosstalk_pixel = zeros(size(d_detector));

    [X_pixel, Y_pixel] = meshgrid(x_pixel, y_pixel);
    for k = 1:length(d_detector)
        shifted_mask = ((X_pixel - d_detector(k))/r0).^2 + ...
                       (Y_pixel/(r0*params.anamorphism)).^2 <= 1;
        crosstalk_pixel(k) = sum(I_pixel(shifted_mask), 'all');
    end
end

%% Dispersed Pixel Integration
% =========================================================================
function [crosstalk_pixel, d_steps, I_pixel, x_pixel, y_pixel] = pixel_integrated_dispersed_analysis(...
    params, lambda_samples, y_pos, r0_values)
% PIXEL_INTEGRATED_DISPERSED_ANALYSIS - Pixel integration for dispersed spectra

    q = params.subsampling_factor;
    nx_pixel = params.fov_pixels;
    x_span = nx_pixel * params.pixel_size;
    nx_sub = nx_pixel * q;
    dx = params.pixel_size / q;
    x_subsampled = ((0:nx_sub-1) + 0.5) * dx - x_span/2;

    margin = max(2 * max(r0_values) * params.anamorphism, 2 * params.pixel_size);
    y_min = min(y_pos) - margin;
    y_max = max(y_pos) + margin;
    ny_pixel = ceil((y_max - y_min) / params.pixel_size);
    y_span = ny_pixel * params.pixel_size;
    y_center = (y_min + y_max) / 2;
    ny_sub = ny_pixel * q;
    dy = params.pixel_size / q;
    y_subsampled = ((0:ny_sub-1) + 0.5) * dy - y_span/2 + y_center;

    [X_subs, Y_subs] = meshgrid(x_subsampled, y_subsampled);
    combined_psf_subsampled = zeros(ny_sub, nx_sub);

    for i = 1:length(lambda_samples)
        R_ellip = sqrt((X_subs/r0_values(i)).^2 + ...
                       ((Y_subs - y_pos(i))/(r0_values(i)*params.anamorphism)).^2);
        PSF = airy_intensity(R_ellip);
        PSF = PSF / sum(PSF(:));
        combined_psf_subsampled = combined_psf_subsampled + PSF;
    end
    combined_psf_subsampled = combined_psf_subsampled / sum(combined_psf_subsampled(:));

    I_pixel = zeros(ny_pixel, nx_pixel);
    for i = 1:nx_pixel
        for j = 1:ny_pixel
            idx_x = (i-1)*q + 1 : i*q;
            idx_y = (j-1)*q + 1 : j*q;
            I_pixel(j,i) = sum(combined_psf_subsampled(idx_y, idx_x), 'all');
        end
    end
    I_pixel = I_pixel / sum(I_pixel(:));

    x_pixel = ((0:nx_pixel-1) + 0.5) * params.pixel_size - x_span/2;
    y_pixel = ((0:ny_pixel-1) + 0.5) * params.pixel_size - y_span/2 + y_center;
    [X_pixel, Y_pixel] = meshgrid(x_pixel, y_pixel);

    d_max = max(30e-6, 1.2 * params.fibre_separation);
    d_steps = linspace(0, d_max, 30);
    d_detector = d_steps * params.M;
    crosstalk_pixel = zeros(size(d_detector));

    for k = 1:length(d_detector)
        shifted_mask = false(size(X_pixel));
        for i = 1:length(lambda_samples)
            mask_i = ((X_pixel - d_detector(k))/r0_values(i)).^2 + ...
                     ((Y_pixel - y_pos(i))/(r0_values(i)*params.anamorphism)).^2 <= 1;
            shifted_mask = shifted_mask | mask_i;
        end
        crosstalk_pixel(k) = sum(I_pixel(shifted_mask), 'all');
    end
end

%% Visualization Functions
% =========================================================================
function fig = visualize_airy_analysis(params, I_continuous, x, crosstalk, d_steps, key_sep, crosstalk_pct, r0)
    fig = figure('Name', 'Airy PSF Crosstalk Analysis', 'Position', [100, 100, 1200, 800]);
    
    % Main PSF visualization
    subplot(2, 2, [1, 3]);
    imagesc(x*1e6, x*1e6, log10(I_continuous + 1e-10));
    colormap hot; colorbar;
    caxis([-7.5, 0]);
    axis image; hold on;
    
    % Add pixel grid
    pixel_size_um = params.pixel_size * 1e6;
    x_range = xlim(); y_range = ylim();
    for x_line = x_range(1):pixel_size_um:x_range(2)
        plot([x_line, x_line], y_range, 'Color', [0.7, 0.7, 0.7, 0.5], 'LineWidth', 0.5);
    end
    for y_line = y_range(1):pixel_size_um:y_range(2)
        plot(x_range, [y_line, y_line], 'Color', [0.7, 0.7, 0.7, 0.5], 'LineWidth', 0.5);
    end
    
    % Airy disk and fibre annotations
    rectangle('Position', [-r0*1e6, -r0*1e6*params.anamorphism, ...
                          2*r0*1e6, 2*r0*1e6*params.anamorphism], ...
              'Curvature', [1, 1], 'EdgeColor', 'b', 'LineWidth', 1.5);
    
    rectangle('Position', [key_sep*1e6 - r0*1e6, -r0*1e6*params.anamorphism, ...
                          2*r0*1e6, 2*r0*1e6*params.anamorphism], ...
              'Curvature', [1, 1], 'EdgeColor', 'g', 'LineWidth', 1.5);
    
    title('Airy PSF and Fibre Overlap');
    xlabel('Spatial Direction [μm]'); ylabel('Spectral Direction [μm]');
    
    % Crosstalk curve
    subplot(2, 2, 2);
    plot(d_steps*1e6, crosstalk*100, 'LineWidth', 2);
    xline(params.fibre_separation*1e6, 'r--', 'LineWidth', 1.5, ...
          'Label', sprintf('Design: %.3f%%', crosstalk_pct));
    xlabel('Fibre Separation [μm]'); ylabel('Crosstalk [%]');
    title('Crosstalk vs Separation'); grid on;
    
    % System parameters
    subplot(2, 2, 4);
    text(0.1, 0.9, sprintf('System Parameters:'), 'FontWeight', 'bold', 'FontSize', 11);
    text(0.1, 0.7, sprintf('Fibre diameter: %.1f μm', params.fibre_diameter*1e6));
    text(0.1, 0.6, sprintf('Magnification: %.2f', params.M));
    text(0.1, 0.5, sprintf('F/# camera: %.2f', params.F2));
    text(0.1, 0.4, sprintf('Anamorphism: %.2f', params.anamorphism));
    text(0.1, 0.3, sprintf('Wavelength: %.1f nm', params.wavelength*1e9));
    text(0.1, 0.2, sprintf('Airy radius: %.1f μm', r0*1e6));
    axis off;
end

function fig = visualize_gaussian_analysis(params, I_gaussian, x, crosstalk, d_steps, key_sep, crosstalk_pct, wx, wy)
    fig = figure('Name', 'Gaussian PSF Crosstalk Analysis', 'Position', [100, 100, 1200, 800]);
    
    % Main PSF visualization
    subplot(2, 2, [1, 3]);
    imagesc(x*1e6, x*1e6, log10(I_gaussian + 1e-10));
    colormap hot; colorbar;
    caxis([-7.5, 0]);
    axis image; hold on;
    
    % Add pixel grid
    pixel_size_um = params.pixel_size * 1e6;
    x_range = xlim(); y_range = ylim();
    for x_line = x_range(1):pixel_size_um:x_range(2)
        plot([x_line, x_line], y_range, 'Color', [0.7, 0.7, 0.7, 0.5], 'LineWidth', 0.5);
    end
    for y_line = y_range(1):pixel_size_um:y_range(2)
        plot(x_range, [y_line, y_line], 'Color', [0.7, 0.7, 0.7, 0.5], 'LineWidth', 0.5);
    end
    
    % Gaussian beam contours
    plot_ellipse(0, 0, wx*1e6, wy*1e6, 'b', '-');
    plot_ellipse(key_sep*1e6, 0, wx*1e6, wy*1e6, 'g', '-');
    
    % Fibre core positions
    fibre_radius = params.fibre_radius_detector;
    rectangle('Position', [-fibre_radius*1e6, -fibre_radius*params.anamorphism*1e6, ...
                          2*fibre_radius*1e6, 2*fibre_radius*params.anamorphism*1e6], ...
              'Curvature', [1, 1], 'EdgeColor', 'b', 'LineWidth', 1, 'LineStyle', '--');
    
    rectangle('Position', [key_sep*1e6 - fibre_radius*1e6, -fibre_radius*params.anamorphism*1e6, ...
                          2*fibre_radius*1e6, 2*fibre_radius*params.anamorphism*1e6], ...
              'Curvature', [1, 1], 'EdgeColor', 'g', 'LineWidth', 1, 'LineStyle', '--');
    
    title('Gaussian PSF and Beam Overlap');
    xlabel('Spatial Direction [μm]'); ylabel('Spectral Direction [μm]');
    
    % Crosstalk curve
    subplot(2, 2, 2);
    plot(d_steps*1e6, crosstalk*100, 'LineWidth', 2, 'Color', [0.8500, 0.3250, 0.0980]);
    xline(params.fibre_separation*1e6, 'r--', 'LineWidth', 1.5, ...
          'Label', sprintf('Design: %.3f%%', crosstalk_pct));
    xlabel('Fibre Separation [μm]'); ylabel('Crosstalk [%]');
    title('Crosstalk vs Separation'); grid on;
    
    % System parameters
    subplot(2, 2, 4);
    text(0.1, 0.9, 'Gaussian Beam Parameters:', 'FontWeight', 'bold', 'FontSize', 11);
    text(0.1, 0.7, sprintf('Waist X: %.1f μm', wx*1e6));
    text(0.1, 0.6, sprintf('Waist Y: %.1f μm', wy*1e6));
    text(0.1, 0.5, sprintf('MFD@1310nm: %.1f μm', params.MFD_1310*1e6));
    text(0.1, 0.4, sprintf('Anamorphism: %.2f', params.anamorphism));
    text(0.1, 0.3, sprintf('Wavelength: %.1f nm', params.wavelength*1e9));
    axis off;

    % Nested helper function for ellipse plotting
    function plot_ellipse(x0, y0, a, b, color, style)
        theta = linspace(0, 2*pi, 100);
        x_ellipse = x0 + a*cos(theta);
        y_ellipse = y0 + b*sin(theta);
        plot(x_ellipse, y_ellipse, 'Color', color, 'LineStyle', style, 'LineWidth', 1.5);
    end
end

function fig = visualize_dispersed_analysis(params, combined_psf, x, y, y_pos, r0_values, lambda_samples, key_sep, crosstalk_pct, d_steps, crosstalk)
    fig = figure('Name', 'Dispersed Spectrum Crosstalk', 'Position', [100, 100, 1200, 800]);
    
    % Main PSF visualization
    subplot(2, 2, [1, 3]);
    imagesc(x*1e6, y*1e6, log10(combined_psf + 1e-10));
    colormap hot; colorbar;
    caxis([-7.5, 0]);
    axis image; hold on;
    
    % Add pixel grid
    pixel_size_um = params.pixel_size * 1e6;
    x_range = xlim(); y_range = ylim();
    for x_line = x_range(1):pixel_size_um:x_range(2)
        plot([x_line, x_line], y_range, 'Color', [0.7, 0.7, 0.7, 0.5], 'LineWidth', 0.5);
    end
    for y_line = y_range(1):pixel_size_um:y_range(2)
        plot(x_range, [y_line, y_line], 'Color', [0.7, 0.7, 0.7, 0.5], 'LineWidth', 0.5);
    end
    
    % Spectral PSF positions - original fibres (blue)
    for i = 1:length(lambda_samples)
        rectangle('Position', [-r0_values(i)*1e6, (y_pos(i)-r0_values(i)*params.anamorphism)*1e6, ...
                              2*r0_values(i)*1e6, 2*r0_values(i)*params.anamorphism*1e6], ...
                  'EdgeColor', 'b', 'LineWidth', 0.5);
    end
    
    % Spectral PSF positions - shifted fibres (green)
    for i = 1:length(lambda_samples)
        rectangle('Position', [key_sep*1e6 - r0_values(i)*1e6, (y_pos(i)-r0_values(i)*params.anamorphism)*1e6, ...
                              2*r0_values(i)*1e6, 2*r0_values(i)*params.anamorphism*1e6], ...
                  'EdgeColor', 'g', 'LineWidth', 0.5);
    end
    
    title('Dispersed Spectrum PSF Distribution');
    xlabel('Spatial Direction [μm]'); ylabel('Spectral Direction [μm]');
    
    % Spectral dispersion plot
    subplot(2, 2, 2);
    plot(lambda_samples*1e9, y_pos*1e6, '-', 'LineWidth', 2);
    hold on;
    plot(lambda_samples(display_idx)*1e9, y_pos(display_idx)*1e6, ...
         'o', 'MarkerSize', 5, 'LineWidth', 1.2);
    hold off;
    xlabel('Wavelength [nm]'); ylabel('Spectral Position [μm]');
    title('Spectral Dispersion'); grid on;
    
    % Crosstalk curve
    subplot(2, 2, 4);
    plot(d_steps*1e6, crosstalk*100, 'LineWidth', 2);
    xline(params.fibre_separation*1e6, 'r--', 'LineWidth', 1.5, ...
          'Label', sprintf('Design: %.3f%%', crosstalk_pct));
    xlabel('Fibre Separation [μm]'); ylabel('Crosstalk [%]');
    title('Crosstalk vs Separation'); grid on;
end

%% Model Comparison Function
% =========================================================================
function comparison = compare_crosstalk_models(crosstalk_results, params)
% COMPARE_CROSSTALK_MODELS - Generate comparative analysis of all models

    fprintf('\n--- Model Comparison ---\n');
    
    % Extract crosstalk values at design separation
    design_sep = params.fibre_separation;
    
    airy_val = crosstalk_results.airy.continuous.key_separation_crosstalk;
    gaussian_val = crosstalk_results.gaussian.key_separation_crosstalk;
    dispersed_val = crosstalk_results.dispersed.continuous.key_separation_crosstalk;
    
    fprintf('Crosstalk at %.1f μm separation:\n', design_sep*1e6);
    fprintf('  Airy model:     %.3f%%\n', airy_val);
    fprintf('  Gaussian model: %.3f%%\n', gaussian_val);
    fprintf('  Dispersed model: %.3f%%\n', dispersed_val);
    
    vals = [airy_val, gaussian_val, dispersed_val];
    mean_val = mean(vals);
    if mean_val > 0
        relative_spread = std(vals) / mean_val;
    else
        relative_spread = NaN;
    end

    comparison = struct(...
        'design_separation', design_sep, ...
        'airy_crosstalk', airy_val, ...
        'gaussian_crosstalk', gaussian_val, ...
        'dispersed_crosstalk', dispersed_val, ...
        'relative_model_spread', relative_spread);

    fprintf('Relative spread between model-specific leakage metrics: %.1f%%\n', ...
        100 * relative_spread);
end

function I = airy_intensity(u)
% AIRY_INTENSITY - Airy intensity with u normalized to the first-zero radius
    z = 1.22 * pi * u;
    I = ones(size(z));
    nz = abs(z) > 1e-12;
    I(nz) = (2 * besselj(1, z(nz)) ./ z(nz)).^2;
end
