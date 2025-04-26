%% Massive MIMO System Simulation
clear all; close all; clc;
rng(42);

%% System Parameters
Nt = 64;          % Number of transmit antennas
Nr = 20;          % Number of receive antennas
N_users = 10;     % Number of users
Nr_per_user = 2;  % Number of receive antennas per user

% URA parameters
Nx = 8;           
Ny = 8;          

N_sc = 128;       % Number of subcarriers
mod_order = 4;    % Modulation order
 
% Channel parameters
K_factor = 3;      % Rician K-factor
path_loss_exp = 4; % Path loss exponent
cell_radius = 500; % Cell radius in meters
min_distance = 35; % Minimum distance from BS

SNR_dB = 0:2:30;
num_symbols = 200; 
num_trials = 2;

% Display system configuration
fprintf('\nSystem Configuration:\n');
fprintf('- %d Tx antennas in %dx%d URA configuration\n', Nt, Nx, Ny);
fprintf('- %d Rx antennas (%d per user, %d users)\n', Nr, Nr_per_user, N_users);
fprintf('- %d OFDM subcarriers \n', N_sc);
fprintf('- Modulation order: %d\n', mod_order);
fprintf('- Cell radius: %d m, Min distance: %d m\n', cell_radius, min_distance);
fprintf('- Number of trials: %d\n', num_trials);

%% Set up URA antenna geometry
[tx_x, tx_y] = meshgrid(0:Ny-1, 0:Nx-1);
tx_positions = [tx_x(:), tx_y(:)];
tx_positions = tx_positions * 0.5;

fc = 28e9; lambda = physconst('LightSpeed')/fc;
array = phased.URA('Size', [Nx, Ny], 'ElementSpacing', lambda/2);
figure; viewArray(array); grid on;
figure; pattern(array, fc, 'PropagationSpeed', physconst('LightSpeed'));
title('Radiation Pattern of URA'); xlabel('Angle (degrees)'); ylabel('Gain (dB)'); grid on;

%% Generate random user positions

user_positions = zeros(N_users, 2);
user_distances = zeros(N_users, 1);
user_angles = zeros(N_users, 1);

for u = 1:N_users
    angle = 2 * pi * rand();
    distance = min_distance + (cell_radius - min_distance) * rand();
    user_positions(u, :) = [distance * cos(angle), distance * sin(angle)];
    user_distances(u) = distance;
    user_angles(u) = angle * 180/pi;
end

figure;
scatter(0, 0, 100, 'filled', 'b', 'DisplayName', 'Base Station'); hold on;
scatter(user_positions(:, 1), user_positions(:, 2), 50, 'filled', 'r', 'DisplayName', 'Users'); grid on;
legend();
title('User and Base Station Positions');
xlabel('x-coordinate (m)'); ylabel('y-coordinate (m)');
axis equal;

%% Rx Power Simulation

viewer = siteviewer('Buildings','1_map.osm', 'Basemap', 'topographic');

bs_lat = 12.843640;
bs_lon = 80.153541;
bs_height = 120;

tx = txsite('Name', 'Base Station', ...
            'Antenna', array, ...
            'Latitude', bs_lat, ...
            'Longitude', bs_lon, ...
            'AntennaHeight', bs_height, ...
            'TransmitterFrequency', fc);

pattern(tx, 'Transparency', 0.6);
mobile_heights = 1.5 * ones(N_users, 1);

rx_sites = rxsite.empty;
for i = 1:N_users
    [rx_lat, rx_lon] = location(tx, user_distances(i), user_angles(i));
    rx_sites(i) = rxsite('Name', sprintf('Mobile %d', i), ...
                'Latitude', rx_lat, ...
                'Longitude', rx_lon, ...
                'AntennaHeight', mobile_heights(i));
    show(rx_sites(i));
end

prop_model = propagationModel('raytracing', ...
                             'MaxNumReflections', 1, ...
                             'MaxNumDiffractions', 0);

rx_power = zeros(N_users, 1);
for i = 1:N_users
    rx_power(i) = sigstrength(rx_sites(i), tx, prop_model);
end

rx_power_table = table((1:N_users)', rx_power, 'VariableNames', {'User', 'RxPower_dBm'});
disp('Received Power at Mobile Stations:');
disp(rx_power_table);
coverage(tx, prop_model);

%% QAM Modulation Configuration

bits_per_symbol = log2(mod_order);
qam_const = qammod(0:mod_order-1, mod_order, 'UnitAveragePower', true);

%% Performance metrics
ber_all = zeros(length(SNR_dB), num_trials);
ser_all = zeros(length(SNR_dB), num_trials);
capacity_all = zeros(length(SNR_dB), num_trials);

performance_metrics = table('Size', [length(SNR_dB), 4], ...
                          'VariableTypes', {'double', 'double', 'double', 'double'}, ...
                          'VariableNames', {'SNR_dB', 'BER', 'SER', 'Capacity_bps_Hz'});
performance_metrics.SNR_dB = SNR_dB';

% main simulation loop
for trial = 1:num_trials
    fprintf('Running trial %d of %d\n', trial, num_trials);
    
    % channel matrix for each subcarrier
    H_freq = zeros(Nr, Nt, N_sc);
    
    % channel generation for each user
    for u = 1:N_users
        user_idx = ((u-1)*Nr_per_user+1):(u*Nr_per_user);
        
        user_distance = user_distances(u);
        user_pos = user_positions(u, :);
        
        % Path loss (large scale fading)
        path_loss = (user_distance/min_distance)^(-path_loss_exp);
        
        for rx_ant = 1:Nr_per_user
            for tx_ant = 1:Nt
                % Position difference between Tx and Rx antenna
                dx = user_pos(1) - tx_positions(tx_ant, 1);
                dy = user_pos(2) - tx_positions(tx_ant, 2);
                distance = sqrt(dx^2 + dy^2 + user_distance^2);
                
                % LOS component
                phase_los = 2*pi*distance;
                h_los = exp(-1j * phase_los);
                
                for sc = 1:N_sc
                    % Rayleigh component
                    h_nlos = (randn() + 1j*randn())/sqrt(2);
                    K = K_factor;
                    H_freq(user_idx(rx_ant), tx_ant, sc) = sqrt(path_loss) * ...
                        sqrt(K/(K+1))*h_los + sqrt(1/(K+1))*h_nlos;
                end
            end
        end
    end
    
    for snr_idx = 1:length(SNR_dB)
        snr_db = SNR_dB(snr_idx);
        snr_linear = 10^(snr_db/10);
        
        total_bits = 0;
        total_errors = 0;
        total_symbols = 0;
        total_symbol_errors = 0;
        
        for sym = 1:num_symbols
            tx_bits = randi([0 1], N_users, N_sc * bits_per_symbol);
            tx_symbols = zeros(N_users, N_sc);
            for u = 1:N_users
                for sc = 1:N_sc
                    bit_idx = (sc-1)*bits_per_symbol + 1;
                    user_bits = tx_bits(u, bit_idx:(bit_idx+bits_per_symbol-1));
                    decimal_value = bi2de(user_bits, 'left-msb');
                    tx_symbols(u, sc) = qam_const(decimal_value + 1);
                end
            end
            
            % channel matrices for each user
            H_users = zeros(N_users, Nt, N_sc);
            for u = 1:N_users
                user_idx = ((u-1)*Nr_per_user+1):(u*Nr_per_user);
                for sc = 1:N_sc
                    H_users(u, :, sc) = mean(H_freq(user_idx, :, sc), 1);
                end
            end
            
            tx_precoded = zeros(Nt, N_sc);
            
            for sc = 1:N_sc
                H_sc = squeeze(H_users(:, :, sc));
                W_zf = H_sc' * inv(H_sc * H_sc'); % zf precoding
                normalization = sqrt(N_users / trace(W_zf * W_zf'));
                W_zf = normalization * W_zf;
                tx_precoded(:, sc) = W_zf * tx_symbols(:, sc);
            end
            
            % awgn
            noise_var = 1 / snr_linear;
            rx_noise = sqrt(noise_var/2) * (randn(Nr, N_sc) + 1j*randn(Nr, N_sc));
            
            % pass through channel
            rx_symbols = zeros(Nr, N_sc);
            for sc = 1:N_sc
                rx_symbols(:, sc) = H_freq(:, :, sc) * tx_precoded(:, sc) + rx_noise(:, sc);
            end
            
            rx_eq_symbols = zeros(N_users, N_sc);
            for u = 1:N_users
                user_idx = ((u-1)*Nr_per_user+1):(u*Nr_per_user);
                user_rx = rx_symbols(user_idx, :);
                
                % simple combining for user with multiple antennas
                for sc = 1:N_sc
                    user_H = H_freq(user_idx, :, sc);
                    
                    % zf equalization
                    H_sc = squeeze(H_users(:, :, sc));
                    W_zf = H_sc' * inv(H_sc * H_sc');
                    normalization = sqrt(N_users / trace(W_zf * W_zf'));
                    W_zf = normalization * W_zf;
                    
                    w_mrc = user_H * W_zf(:, u);
                    w_mrc = w_mrc / norm(w_mrc);
                    rx_eq_symbols(u, sc) = (w_mrc' * user_rx(:, sc));
                end
            end
            
            rx_bits = zeros(N_users, N_sc * bits_per_symbol);
            for u = 1:N_users
                for sc = 1:N_sc
                    rx_sym = rx_eq_symbols(u, sc);
                    
                    [~, dec_symbol] = min(abs(rx_sym - qam_const));
                    dec_symbol = dec_symbol - 1;
                    
                    bit_idx = (sc-1)*bits_per_symbol + 1;
                    demod_bits = de2bi(dec_symbol, bits_per_symbol, 'left-msb');
                    rx_bits(u, bit_idx:(bit_idx+bits_per_symbol-1)) = demod_bits;
                    
                    orig_dec = bi2de(tx_bits(u, bit_idx:(bit_idx+bits_per_symbol-1)), 'left-msb');
                    
                    if dec_symbol ~= orig_dec
                        total_symbol_errors = total_symbol_errors + 1;
                    end
                    total_symbols = total_symbols + 1;
                end
                
                bit_errors = sum(tx_bits(u, :) ~= rx_bits(u, :));
                total_errors = total_errors + bit_errors;
                total_bits = total_bits + length(tx_bits(u, :));
            end
        end
        
        ber = total_errors / total_bits;
        ser = total_symbol_errors / total_symbols;
        ber_all(snr_idx, trial) = ber;
        ser_all(snr_idx, trial) = ser;
        
        capacity = 0;
        for sc = 1:N_sc
            H = H_freq(:, :, sc);
            capacity = capacity + log2(det(eye(Nr) + snr_linear/Nt * (H * H')));
        end
        capacity = real(capacity) / N_sc;
        capacity_all(snr_idx, trial) = capacity;
        
        performance_metrics.BER(snr_idx) = ber;
        performance_metrics.SER(snr_idx) = ser;
        performance_metrics.Capacity_bps_Hz(snr_idx) = capacity;
    end
end

disp('Performance Metrics:');
disp(performance_metrics);

%% Average results over trials
ber_avg = mean(ber_all, 2);
ser_avg = mean(ser_all, 2);
capacity_avg = mean(capacity_all, 2);

%% Plot
figure;
plot(SNR_dB, capacity_avg, 'd-', 'LineWidth', 2);
grid on; xlim ([0 20]);
xlabel('SNR (dB)');
ylabel('Capacity (bps/Hz)');
title('Capacity of Massive MIMO System');