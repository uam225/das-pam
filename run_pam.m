% --- Load the recorded data ---
load('/Users/umairarshad/MRes/kWave_simulation/sim_data/simulation_ppw10.0_2025-03-29_15-55-25.mat');  % This should load `sensor_data`

% --- Define simulation and array parameters ---
c = 1500;               % Speed of sound [m/s] → 1.5 mm/us
Nx = 667;
Ny = 667;
dx = 100 / Nx;          % Grid spacing [mm]
dt = 0.5 * (dx / c);    % Time step [us] (matches simulation CFL condition)
fprintf('Time step (dt): %.6f µs\n', dt);

% --- Construct grid vectors (covering full 100 mm x 100 mm domain) ---
grid_x = (0:Nx-1) * dx;  % 0 to 99.5 mm
grid_y = (0:Ny-1) * dx;

% --- Reconstruct sensor positions based on your simulation setup ---
num_receivers = 64;
x_sensor = linspace(10, 90, num_receivers);  % 10 mm to 90 mm
y_sensor = ones(1, num_receivers) * 90;      % All at y = 90 mm
sensor_positions = [x_sensor', y_sensor'];


% --- Call the PCM beamformer ---
tic
[pcm_map, x_scan, y_scan] = generate_pam(sensor_data, sensor_positions, grid_x, grid_y, c, dt);
disp('Beamformer completed')
toc
%%
% --- Visualize the result ---
figure;
imagesc(x_scan, y_scan, pcm_map);
axis image;
xlabel('X [mm]');
ylabel('Y [mm]');
title('Cavitation Map');
colorbar;
set(gca, 'YDir', 'normal');  % Flip Y axis for correct orientation

