function [pcm_map, x_scan, y_scan] = generate_pam(sensor_data, sensor_positions, grid_x, grid_y, c, dt)
% generate_pcm - Passive Cavitation Mapping using delay-and-sum beamforming
%
% INPUTS:
%   sensor_data      : [nChannels x nTime] matrix of time-series data
%   sensor_positions : [nChannels x 2] positions of each sensor in mm [x, y]
%   grid_x           : vector of x positions for the image grid [1 x Nx] in mm
%   grid_y           : vector of y positions for the image grid [1 x Ny] in mm
%   c                : speed of sound in mm/us
%   dt               : time step between samples in us
%
% OUTPUT:
%   pcm_map          : [Ny x Nx] matrix of source power estimates

% Initialization
[nChannels, nSamples] = size(sensor_data);
Nx = length(grid_x);
Ny = length(grid_y);
dx = 100 / Nx;
dy = dx;
%pcm_map = zeros(Ny, Nx);                % Final power map
time_vector = (0:nSamples-1) * dt;      % Time axis

step = 2;  % Try 2 first. Increase to 3 or more for faster tests

x_scan = 0:step*dx:(Nx-1)*dx;
y_scan = 0:step*dy:(Ny-1)*dy;

Nx_scan = length(x_scan);
Ny_scan = length(y_scan);

pcm_map = zeros(Nx_scan, Ny_scan);
% Loop over each pixel in the image grid
for ix = 1:Nx_scan
    for iy = 1:Ny_scan
        % Get pixel position
        px = x_scan(ix);
        py = y_scan(iy);
        delayed_signals = zeros(nChannels, nSamples);  % Preallocate % Try flipping nSamples and nChannels to fix orientation difference

        for ch = 1:nChannels
            % Sensor position
            sx = sensor_positions(ch,1);
            sy = sensor_positions(ch,2);

            % Compute distance and time delay
            dist = sqrt((px - sx)^2 + (py - sy)^2);   % mm
            delay = dist / c;                         % us
           
            % Interpolate signal at delayed time
            delayed_signals(ch, :) = interp1(time_vector, ...
                                            sensor_data(ch,:), ...
                                            time_vector + delay, ...
                                            'linear', 0);
        end

        % Delay-and-sum across all channels
        summed_signal = sum(delayed_signals, 1);

        % Estimate source power (squared amplitude over time)
        pcm_map(iy, ix) = sum(summed_signal.^2);
    end
end

end