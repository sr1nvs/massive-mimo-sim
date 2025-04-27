%% Average results over trials
ber_avg = mean(ber_all, 2);
ser_avg = mean(ser_all, 2);
capacity_avg = mean(capacity_all, 2);

%% Plot
figure;
plot(SNR_dB, capacity_avg, 'd-', 'LineWidth', 2);
grid on; xlim([0 20]);
xlabel('SNR (dB)');
ylabel('Capacity (bps/Hz)');
title('Capacity of Massive MIMO System with 5% Channel Estimation Error');

figure;
semilogy(SNR_dB, ber_avg, 'o-', 'LineWidth', 2);
grid on; xlim([0 20]); ylim([1e-6 1]);
xlabel('SNR (dB)');
ylabel('BER');
title('BER Performance with 5% Channel Estimation Error');

figure;
semilogy(SNR_dB, ser_avg, 's-', 'LineWidth', 2);
grid on; xlim([0 20]); ylim([1e-6 1]);
xlabel('SNR (dB)');
ylabel('SER');
title('SER Performance with 5% Channel Estimation Error');