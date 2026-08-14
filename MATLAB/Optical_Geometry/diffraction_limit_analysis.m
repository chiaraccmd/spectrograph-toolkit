function [transition_data, performance_metrics] = diffraction_limit_analysis(varargin)
% DIFFRACTION_LIMIT_ANALYSIS - Determine spectrograph performance transition points
%
% Analyzes when a spectrograph transitions from geometric (slit-limited) to
% diffraction-limited performance. Identifies operating regimes and
% performance boundaries for optical design.
%
% Physical principles:
% - Geometric resolving power: R_geom = (m·G·λ·F·W) / s1
% - Diffraction-limited resolving power: R_diff = (m·G·W) / 1.22
% - Transition occurs when R_geom = R_diff
%
% Inputs (parameter/value pairs):
%   'wavelength_range' - [min,max] wavelength [m] (default: [0.3e-6, 1.7e-6])
%   'grating_density'  - Grating lines per meter (default: 650e3)
%   'beam_size'        - Beam size at grating [m] (default: 14.8e-3)
%   'f_number'         - Collimator focal ratio (default: 3.57)
%   'slit_width'       - Slit width [m] (default: 4.14e-6)
%   'diffraction_order'- Diffraction order (default: 1)
%   'n_points'         - Number of wavelength points (default: 1000)
%
% Outputs:
%   transition_data    - Structure with transition wavelength and performance data
%   performance_metrics- Structure with key performance indicators
%
% Example:
%   % Basic analysis with default parameters
%   [transition, metrics] = diffraction_limit_analysis();
%
%   % Custom design study
%   [transition, metrics] = diffraction_limit_analysis(...
%       'grating_density', 800e3, 'beam_size', 20e-3, 'f_number', 4.0);

%% Input parsing and parameter setup
% =========================================================================
p = inputParser;

addParameter(p, 'wavelength_range', [0.3e-6, 1.7e-6], @(x) numel(x)==2 && all(x>0));
addParameter(p, 'grating_density', 650e3, @isnumeric);      % Grating lines/m
addParameter(p, 'beam_size', 14.8e-3, @isnumeric);          % Beam size [m]
addParameter(p, 'f_number', 3.57, @isnumeric);              % Collimator f-number
addParameter(p, 'slit_width', 4.14e-6, @isnumeric);         % Slit width [m]
addParameter(p, 'diffraction_order', 1, @isnumeric);        % Diffraction order
addParameter(p, 'n_points', 1000, @isnumeric);              % Wavelength points

parse(p, varargin{:});
params = p.Results;

% Create wavelength array
lambda = linspace(params.wavelength_range(1), params.wavelength_range(2), params.n_points);

fprintf('=== Spectrograph Diffraction Limit Analysis ===\n');
fprintf('Grating: %.0f l/mm, Beam: %.1f mm, F/#: %.2f, Slit: %.1f μm\n', ...
    params.grating_density*1e-3, params.beam_size*1e3, params.f_number, params.slit_width*1e6);

%% Resolving Power Calculations
% =========================================================================
% Geometric (slit-limited) resolving power
% R_geom = (m * G * λ * F * W) / s1
R_geometric = (params.diffraction_order * params.grating_density .* lambda * ...
               params.f_number * params.beam_size) ./ params.slit_width;

% Diffraction-limited resolving power (Rayleigh criterion)
% R_diff = (m * G * W) / 1.22
R_diffraction = params.diffraction_order * params.grating_density * params.beam_size / 1.22;

fprintf('Diffraction-limited resolving power: R = %.0f\n', R_diffraction);

%% Transition Point Analysis
% =========================================================================
% Find where geometric limit equals diffraction limit
transition_idx = find(R_geometric >= R_diffraction, 1);

if ~isempty(transition_idx)
    lambda_transition = lambda(transition_idx);
    R_transition = R_diffraction;
    
    fprintf('Transition wavelength: %.2f μm\n', lambda_transition*1e6);
    fprintf('At λ < %.2f μm: system is slit-limited\n', lambda_transition*1e6);
    fprintf('At λ > %.2f μm: system is diffraction-limited\n', lambda_transition*1e6);
else
    lambda_transition = NaN;
    R_transition = NaN;
    fprintf('System remains slit-limited across entire wavelength range\n');
end

%% Performance Region Analysis
% =========================================================================
% Define performance regions
slit_limited_mask = R_geometric < R_diffraction;
diffraction_limited_mask = R_geometric >= R_diffraction;

% Calculate performance metrics
if any(slit_limited_mask)
    max_geometric_R = max(R_geometric(slit_limited_mask));
else
    max_geometric_R = NaN;
end

if ~isnan(max_geometric_R)
    fprintf('Maximum geometric (slit-limited) R: %.0f\n', max_geometric_R);
end

%% Visualization
% =========================================================================
fig = create_performance_plot(lambda, R_geometric, R_diffraction, ...
    lambda_transition, R_transition, params);

%% Performance Optimization Analysis
% =========================================================================
% Calculate how changing parameters affects transition point
optimization_data = analyze_parameter_sensitivity(params);

%% Prepare Output Structures
% =========================================================================
transition_data = struct(...
    'transition_wavelength', lambda_transition, ...
    'transition_resolving_power', R_transition, ...
    'wavelength_array', lambda, ...
    'geometric_resolving_power', R_geometric, ...
    'diffraction_resolving_power', R_diffraction * ones(size(lambda)), ...
    'performance_regions', struct(...
        'slit_limited', slit_limited_mask, ...
        'diffraction_limited', diffraction_limited_mask));

performance_metrics = struct(...
    'max_geometric_R', max_geometric_R, ...
    'diffraction_limit_R', R_diffraction, ...
    'parameter_sensitivity', optimization_data, ...
    'design_parameters', params);

fprintf('\n=== Analysis Complete ===\n');

end

%% Visualization Function
% =========================================================================
function fig = create_performance_plot(lambda, R_geom, R_diff, lambda_trans, R_trans, params)
% CREATE_PERFORMANCE_PLOT - Generate comprehensive performance visualization

    fig = figure('Name', 'Spectrograph Performance Transition', ...
        'Position', [100, 100, 1200, 800]);
    
    % Main resolving power plot
    subplot(2, 2, [1, 3]);
    hold on;
    
    % Plot geometric resolving power
    h1 = plot(lambda * 1e6, R_geom, 'b-', 'LineWidth', 2.5, ...
        'DisplayName', 'Geometric Limit (Slit-Limited)');
    
    % Plot diffraction limit
    h2 = plot([min(lambda), max(lambda)] * 1e6, [R_diff, R_diff], 'r--', ...
        'LineWidth', 2, 'DisplayName', 'Diffraction Limit');
    
    % Mark transition point if it exists
    if ~isnan(lambda_trans)
        h3 = plot(lambda_trans * 1e6, R_trans, 'ko', ...
            'MarkerSize', 8, 'MarkerFaceColor', 'k', ...
            'DisplayName', 'Transition Point');
        
        % Add vertical line at transition
        plot([lambda_trans, lambda_trans] * 1e6, [0, R_trans], 'k:', ...
            'LineWidth', 1, 'HandleVisibility', 'off');
        
        % Shade performance regions
        y_limits = ylim;
        fill([min(lambda), lambda_trans, lambda_trans, min(lambda)] * 1e6, ...
             [y_limits(1), y_limits(1), y_limits(2), y_limits(2)], ...
             [0.8, 0.9, 1.0], 'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
             'DisplayName', 'Slit-Limited Region');
        
        fill([lambda_trans, max(lambda), max(lambda), lambda_trans] * 1e6, ...
             [y_limits(1), y_limits(1), y_limits(2), y_limits(2)], ...
             [1.0, 0.9, 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
             'DisplayName', 'Diffraction-Limited Region');
    end
    
    xlabel('Wavelength [μm]', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Resolving Power R', 'FontSize', 12, 'FontWeight', 'bold');
    title('Spectrograph Performance: Geometric vs Diffraction Limit', ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    legend('show', 'Location', 'northwest');
    grid on;
    set(gca, 'FontSize', 11);
    
    % Add transition annotation
    if ~isnan(lambda_trans)
        annotation('textbox', [0.15, 0.7, 0.2, 0.1], 'String', ...
            sprintf('Transition:\nλ = %.2f μm\nR = %.0f', ...
            lambda_trans*1e6, R_trans), ...
            'BackgroundColor', 'white', 'FontSize', 10, ...
            'EdgeColor', 'black', 'FitBoxToText', 'on');
    end
    
    %% Parameter sensitivity subplot
    subplot(2, 2, 2);
    
    % Calculate how transition changes with slit width
    slit_widths = linspace(params.slit_width * 0.5, params.slit_width * 2, 50);
    transition_wavelengths = zeros(size(slit_widths));
    
    for i = 1:length(slit_widths)
        R_geom_temp = (params.diffraction_order * params.grating_density .* lambda * ...
                      params.f_number * params.beam_size) ./ slit_widths(i);
        idx = find(R_geom_temp >= R_diff, 1);
        if ~isempty(idx)
            transition_wavelengths(i) = lambda(idx) * 1e6;
        else
            transition_wavelengths(i) = NaN;
        end
    end
    
    plot(slit_widths * 1e6, transition_wavelengths, 'LineWidth', 2, ...
         'Color', [0.2, 0.6, 0.2]);
    xlabel('Slit Width [μm]', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Transition Wavelength [μm]', 'FontSize', 11, 'FontWeight', 'bold');
    title('Slit Width vs Transition Point', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    
    % Mark current design point
    hold on;
    plot(params.slit_width * 1e6, lambda_trans * 1e6, 'ro', ...
         'MarkerSize', 8, 'MarkerFaceColor', 'r', ...
         'DisplayName', 'Current Design');
    
    %% System parameters display
    subplot(2, 2, 4);
    
    param_text = {
        sprintf('System Parameters:')
        sprintf('Grating density: %.0f l/mm', params.grating_density*1e-3)
        sprintf('Beam size: %.1f mm', params.beam_size*1e3)
        sprintf('F-number: %.2f', params.f_number)
        sprintf('Slit width: %.2f μm', params.slit_width*1e6)
        sprintf('Diffraction order: %d', params.diffraction_order)
        ''
        sprintf('Performance:')
        sprintf('Diffraction limit: R = %.0f', R_diff)
    };
    
    if ~isnan(lambda_trans)
        param_text{end+1} = sprintf('Transition: %.2f μm', lambda_trans*1e6);
        geom_values = R_geom(R_geom < R_diff);
        if isempty(geom_values)
            max_geom_R = NaN;
        else
            max_geom_R = max(geom_values);
        end
        param_text{end+1} = sprintf('Max geometric R: %.0f', max_geom_R);
    else
        param_text{end+1} = 'Always slit-limited';
        param_text{end+1} = sprintf('Max R: %.0f', max(R_geom));
    end
    
    text(0.05, 0.95, param_text, 'VerticalAlignment', 'top', ...
         'FontSize', 10, 'FontName', 'FixedWidth', ...
         'BackgroundColor', [0.95, 0.95, 0.95], ...
         'EdgeColor', 'black', 'Margin', 10);
    axis off;
end

%% Parameter Sensitivity Analysis
% =========================================================================
function sensitivity_data = analyze_parameter_sensitivity(params)
% ANALYZE_PARAMETER_SENSITIVITY - Analyze how design changes affect transition

    lambda_test = linspace(params.wavelength_range(1), params.wavelength_range(2), 500);

    % Test parameter variations
    param_variations = struct(...
        'slit_width', linspace(params.slit_width * 0.3, params.slit_width * 3, 20), ...
        'f_number', linspace(params.f_number * 0.5, params.f_number * 2, 20), ...
        'grating_density', linspace(params.grating_density * 0.5, params.grating_density * 2, 20), ...
        'beam_size', linspace(params.beam_size * 0.5, params.beam_size * 2, 20));
    
    sensitivity_data = struct();
    param_names = fieldnames(param_variations);
    
    for i = 1:length(param_names)
        param_name = param_names{i};
        param_values = param_variations.(param_name);
        transition_points = zeros(size(param_values));
        
        for j = 1:length(param_values)
            % Create temporary parameters
            temp_params = params;
            temp_params.(param_name) = param_values(j);
            
            % Calculate geometric resolving power
            R_geom_temp = (temp_params.diffraction_order * temp_params.grating_density .* lambda_test * ...
                          temp_params.f_number * temp_params.beam_size) ./ temp_params.slit_width;
            
            % Diffraction limit must be recomputed for the varied design
            R_diff_temp = temp_params.diffraction_order * temp_params.grating_density * ...
                          temp_params.beam_size / 1.22;
            
            % Find transition point
            idx = find(R_geom_temp >= R_diff_temp, 1);
            if ~isempty(idx)
                transition_points(j) = lambda_test(idx) * 1e6;
            else
                transition_points(j) = NaN;
            end
        end
        
        % Calculate sensitivity
        valid_indices = ~isnan(transition_points);
        if sum(valid_indices) > 1
            valid_transitions = transition_points(valid_indices);
            valid_params = param_values(valid_indices);
            
            % Calculate mean sensitivity across the range
            transition_diff = diff(valid_transitions);
            param_diff = diff(valid_params);
            sensitivity_val = mean(transition_diff ./ param_diff);
        else
            sensitivity_val = NaN;
        end
        
        sensitivity_data.(param_name) = struct(...
            'values', param_values, ...
            'transition_wavelengths', transition_points, ...
            'sensitivity', sensitivity_val);
    end
    
    % Display sensitivity summary
    fprintf('\n--- Parameter Sensitivity Analysis ---\n');
    for i = 1:length(param_names)
        if ~isnan(sensitivity_data.(param_names{i}).sensitivity)
            fprintf('%s sensitivity: %.3f μm per unit change\n', ...
                param_names{i}, sensitivity_data.(param_names{i}).sensitivity);
        end
    end
end
