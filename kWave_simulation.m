clear; clc;

% ---------------------
% DEFINED SETTINGS
% ---------------------
ppw = 15;                        
f_c = 1e6;                       % Center frequency (Hz)
N_c = 5;                        % Number of cycles
domain_length = 100e-3;          % Domain size in meters (100 mm)
c_0 = 1500;                      % Speed of sound (m/s)
density = 1040;                  % Medium density (kg/m³)
CFL = 0.5;                       
num_receivers = 64;              % Number of sensors
sensor_offset = 10e-3;           % Distance from boundary in meters

% ---------------------
% DERIVED PARAMETERS
% ---------------------
lambda = c_0 / f_c;             % Wavelength
dx = lambda / ppw;              % Grid spacing
Nx = round(domain_length / dx); % Grid size in x
dx = domain_length / Nx;        % Adjust dx
Ny = Nx; dy = dx;

dt = CFL * dx / c_0;            % Time step
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% Total simulation time
T_pulse = N_c / f_c;
T_travel = (domain_length - sensor_offset) / c_0;
simulation_time = T_pulse + T_travel;
Nt = round(simulation_time / dt);
kgrid.t_array = (0:Nt-1) * dt;

% ---------------------
% MEDIUM
% ---------------------
medium.sound_speed = c_0;
medium.density = density;

% ---------------------
% SOURCE
% ---------------------
source.p_mask = zeros(Nx, Ny);
source.p_mask(round(Nx/2), round(Ny/2)) = 1;
%source.p_mask(round(Nx/2.2), round(Ny/2)) = 1;

% Sinusoidal burst (truncated after N_c cycles)
source_signal = sin(2 * pi * f_c * kgrid.t_array);
source_signal((1:length(kgrid.t_array)) > (N_c / f_c) / dt) = 0;
source.p = source_signal;

% ---------------------
% SENSORS
% ---------------------
sensor_x = linspace(sensor_offset, domain_length - sensor_offset, num_receivers);
sensor_y = domain_length - sensor_offset;

sensor_x_grid = round(sensor_x / dx);
sensor_y_grid = round(sensor_y / dy);

sensor.mask = zeros(Nx, Ny);
sensor.mask(sub2ind([Nx, Ny], sensor_x_grid, repmat(sensor_y_grid, size(sensor_x_grid)))) = 1;

% ---------------------
% SIMULATION
% ---------------------
fprintf('Running simulation with PPW = %.1f (Grid: %dx%d)...\n', ppw, Nx, Ny);
tic;
input_args = {'PMLSize', 20, 'PMLAlpha', 2, 'PlotSim', true}; % set PlotSim false if batch running
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
runtime = toc;
fprintf('Simulation completed in %.2f seconds.\n', runtime);

% ---------------------
% METRICS
% ---------------------
peak_pos = max(sensor_data(:));
peak_neg = min(sensor_data(:));
rms_val = sqrt(mean(sensor_data(:).^2));
energy = sum(sum(sensor_data.^2, 1)) * dt;

fprintf('Peak Positive: %.4f | Peak Negative: %.4f | RMS: %.4f | Energy: %.4e\n', ...
    peak_pos, peak_neg, rms_val, energy);
%% analysis
imagesc(sensor_data);
colorbar;
xlabel('time point');
ylabel('receiver number');
% ---------------------
% SAVE DATA
% ---------------------
% Define source location in mm
source_x = round((Nx/2) * dx * 1e3); 
source_y = round((Ny/4) * dy * 1e3);
%source2_x = round((Nx/2.2) * dx * 1e3); 

source_str = sprintf('%dmm_%dmm', source_x, source_y);

% Define timestamp in specified format
timestamp = datestr(datetime('now'), 'mm-dd-HH-MM');

% Generate filename
filename = sprintf('%s_ppw%.0f_%s.mat', source_str, ppw, timestamp);
%save(filename, 'sensor_data', 'kgrid', 'dx', 'dt', 'ppw');
fprintf('Saved data to: %s\n', filename);
