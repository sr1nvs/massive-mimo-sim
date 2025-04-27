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
    
    % channel matrix
    H_freq = zeros(Nr, Nt, N_sc);
    
    % channel generation for each user
    for u = 1:N_users
        user_idx = ((u-1)*Nr_per_user+1):(u*Nr_per_user);
        
        user_distance = user_distances(u);
        user_pos = user_positions(u, :);
        
        % path loss (large scale fading)
        path_loss = (user_distance/min_distance)^(-path_loss_exp);
        
        for rx_ant = 1:Nr_per_user
            for tx_ant = 1:Nt

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
    

    H_est = zeros(size(H_freq));
    for sc = 1:N_sc
        error_matrix = (randn(Nr, Nt) + 1j*randn(Nr, Nt))/sqrt(2);
        H_est(:,:,sc) = H_freq(:,:,sc) + (est_error_percent/100) * abs(H_freq(:,:,sc)) .* error_matrix;
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
            

            H_users_est = zeros(N_users, Nt, N_sc);
            for u = 1:N_users
                user_idx = ((u-1)*Nr_per_user+1):(u*Nr_per_user);
                for sc = 1:N_sc
                    H_users_est(u, :, sc) = mean(H_est(user_idx, :, sc), 1);
                end
            end
            
            tx_precoded = zeros(Nt, N_sc);
            
            for sc = 1:N_sc
                H_sc = squeeze(H_users_est(:, :, sc));
                W_zf = H_sc' * inv(H_sc * H_sc'); % zf precoding with estimated channel
                normalization = sqrt(N_users / trace(W_zf * W_zf'));
                W_zf = normalization * W_zf;
                tx_precoded(:, sc) = W_zf * tx_symbols(:, sc);
            end
            
            % awgn
            noise_var = 1 / snr_linear;
            rx_noise = sqrt(noise_var/2) * (randn(Nr, N_sc) + 1j*randn(Nr, N_sc));
            
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
                    user_H_est = H_est(user_idx, :, sc);
                    
                    % zf equalization with estimated channel
                    H_sc = squeeze(H_users_est(:, :, sc));
                    W_zf = H_sc' * inv(H_sc * H_sc');
                    normalization = sqrt(N_users / trace(W_zf * W_zf'));
                    W_zf = normalization * W_zf;
                    
                    w_mrc = user_H_est * W_zf(:, u);
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
            H = H_est(:, :, sc);
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