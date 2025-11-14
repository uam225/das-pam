clear; clc;

% Define parameter ranges for grid convergence study (Appendix B)
Nx_values = 100:100:1000;      % Spatial grid resolution: 1.5 to 15 PPW
CFL_values = [0.3, 0.5];       % Courant-Friedrichs-Lewy numbers evaluated
N_c_values = 5;                % Five-cycle sinusoidal tone burst

domain_length = 100e-3;        % 100 mm square domain
density_default = 1040;        % Soft-tissue ambient density (kg/m^3)
c_0 = 1500;                    % Speed of sound in medium (m/s)
f_c = 1e6;                     % 1 MHz central frequency (wavelength = 1.5 mm)
num_receivers = 64;            % 64-element linear sensor array
central_sensor = round(num_receivers / 2);

% Points per wavelength corresponding to each Nx
ppw_values = (c_0 / f_c) ./ (domain_length ./ Nx_values);

% Initialise storage arrays
num_nx = length(Nx_values);
num_cfl = length(CFL_values);
num_nc = length(N_c_values);

max_sensor_values = zeros(num_nx, num_cfl, num_nc);
min_sensor_values = zeros(num_nx, num_cfl, num_nc);
mean_sensor_values = zeros(num_nx, num_cfl, num_nc);
total_energy = zeros(num_nx, num_cfl, num_nc);
runtime = zeros(num_nx, num_cfl, num_nc);
correlation = NaN(num_nx, num_cfl, num_nc);

% Store central-sensor waveforms and time vectors for reference comparison
waveforms = cell(num_nx, num_cfl, num_nc);
time_axes = cell(num_nx, num_cfl, num_nc);

% Loop over spatial and temporal discretisation parameters
for i = 1:num_nx
    Nx = Nx_values(i);
    dx = domain_length / Nx;
    Ny = Nx;
    dy = dx;
    kgrid = kWaveGrid(Nx, dx, Ny, dy);

    for j = 1:num_cfl
        CFL = CFL_values(j);
        dt = CFL * dx / c_0;

        for k = 1:num_nc
            N_c = N_c_values(k);

            % Compute required total propagation time dynamically
            T_pulse = N_c / f_c;
            T_travel = 90e-3 / c_0;
            simulation_time = T_pulse + T_travel;
            Nt = round(simulation_time / dt);
            t_array = (0:Nt-1) * dt;

            tic;

            % Medium properties
            medium.sound_speed = c_0;
            medium.density = density_default;

            % Point source at domain centre (50 mm, 50 mm)
            source.p_mask = zeros(Nx, Ny);
            source.p_mask(round(Nx/2), round(Ny/2)) = 1;

            source_signal = sin(2 * pi * f_c * t_array);
            source_signal((1:length(t_array)) > (N_c / f_c) / dt) = 0;
            source.p = source_signal;
            kgrid.t_array = t_array;

            % 64-element linear receiver array along y = 90 mm, x in [10, 90] mm
            sensor.mask = zeros(Nx, Ny);
            mm_to_grid = @(mm) round(mm / dx);
            sensor_x_positions_mm = linspace(10e-3, domain_length - 10e-3, num_receivers);
            sensor_x_positions = mm_to_grid(sensor_x_positions_mm);
            sensor_y_position_mm = domain_length - 10e-3;
            sensor_y_position = mm_to_grid(sensor_y_position_mm);
            sensor_x_positions = max(1, min(sensor_x_positions, Nx));
            sensor_y_position = max(1, min(sensor_y_position, Ny));
            sensor.mask(sub2ind([Nx, Ny], sensor_x_positions, repmat(sensor_y_position, size(sensor_x_positions)))) = 1;

            % Run k-Wave simulation with 20-point PML
            input_args = {'PMLSize', 20, 'PMLAlpha', 2, 'PlotSim', false, 'RecordMovie', false};
            sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
            runtime(i, j, k) = toc;

            % Compute amplitude and time-integrated energy metrics
            max_sensor_values(i, j, k) = max(sensor_data(:));
            min_sensor_values(i, j, k) = min(sensor_data(:));
            mean_sensor_values(i, j, k) = mean(sensor_data(:));
            total_energy(i, j, k) = sum(sensor_data(:).^2) * dt;

            % Save central sensor waveform for post-sweep correlation against Nx = 1000 reference
            waveforms{i, j, k} = double(sensor_data(central_sensor, :));
            time_axes{i, j, k} = t_array;

            fprintf('Completed Nx=%d (PPW=%.1f), CFL=%.1f, N_c=%d in %.2f sec\n', ...
                Nx, ppw_values(i), CFL, N_c, runtime(i, j, k));
        end
    end
end

% Compute Pearson waveform similarity against the high-resolution reference (Nx = 1000)
for j = 1:num_cfl
    for k = 1:num_nc
        ref_wave = waveforms{num_nx, j, k};
        ref_time = time_axes{num_nx, j, k};
        for i = 1:num_nx
            cand_wave = interp1(time_axes{i, j, k}, waveforms{i, j, k}, ref_time, 'linear', 0);
            correlation(i, j, k) = corr(cand_wave(:), ref_wave(:));
        end
    end
end

%% Plot Convergence Metrics Across Grid Resolution (Nx and PPW)
figure;

% (a) Peak Positive Pressure
subplot(2, 2, 1);
plot(Nx_values, squeeze(max_sensor_values(:, 1, 1)), '-o', 'LineWidth', 1.8);
hold on;
plot(Nx_values, squeeze(max_sensor_values(:, 2, 1)), '-s', 'LineWidth', 1.8);
hold off;
xlabel('Grid Resolution (N_x)');
ylabel('Peak Positive Pressure (Pa)');
title('(a) Peak Positive Pressure Convergence');
legend('CFL = 0.3', 'CFL = 0.5', 'Location', 'northeast');
grid on;

% (b) Time-Integrated Total Energy
subplot(2, 2, 2);
plot(Nx_values, squeeze(total_energy(:, 1, 1)), '-o', 'LineWidth', 1.8);
hold on;
plot(Nx_values, squeeze(total_energy(:, 2, 1)), '-s', 'LineWidth', 1.8);
hold off;
xlabel('Grid Resolution (N_x)');
ylabel('Total Energy (Pa^2 \cdot s)');
title('(b) Total Energy Convergence');
legend('CFL = 0.3', 'CFL = 0.5', 'Location', 'northeast');
grid on;

% (c) Waveform Similarity vs High-Resolution Reference (Nx = 1000)
subplot(2, 2, 3);
plot(Nx_values, squeeze(correlation(:, 1, 1)), '-o', 'LineWidth', 1.8);
hold on;
plot(Nx_values, squeeze(correlation(:, 2, 1)), '-s', 'LineWidth', 1.8);
hold off;
xlabel('Grid Resolution (N_x)');
ylabel('Pearson Correlation (vs N_x = 1000)');
title('(c) Waveform Similarity');
legend('CFL = 0.3', 'CFL = 0.5', 'Location', 'southeast');
grid on;

% (d) Computational Runtime
subplot(2, 2, 4);
plot(Nx_values, squeeze(runtime(:, 1, 1)), '-o', 'LineWidth', 1.8);
hold on;
plot(Nx_values, squeeze(runtime(:, 2, 1)), '-s', 'LineWidth', 1.8);
hold off;
xlabel('Grid Resolution (N_x)');
ylabel('Runtime (s)');
title('(d) Computational Cost');
legend('CFL = 0.3', 'CFL = 0.5', 'Location', 'northwest');
grid on;

%% Display Convergence Summary Table (CFL = 0.5)
convergence_table = table( ...
    Nx_values', ...
    ppw_values', ...
    squeeze(max_sensor_values(:, 2, 1)), ...
    squeeze(total_energy(:, 2, 1)), ...
    squeeze(correlation(:, 2, 1)), ...
    squeeze(runtime(:, 2, 1)), ...
    'VariableNames', {'Nx', 'PPW', 'Max_Pressure', 'Total_Energy', 'Waveform_Correlation', 'Runtime_s'} ...
);
disp(convergence_table);
