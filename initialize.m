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