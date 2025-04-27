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