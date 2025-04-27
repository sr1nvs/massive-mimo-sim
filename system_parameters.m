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
K_factor = 12;      % Rician K-factor
path_loss_exp = 4; % Path loss exponent
cell_radius = 500; % Cell radius in meters
min_distance = 35; % Minimum distance from BS

SNR_dB = 0:2:20;
num_symbols = 500; 
num_trials = 1;

% system configuration
fprintf('\nSystem Configuration:\n');
fprintf('- %d Tx antennas in %dx%d URA configuration\n', Nt, Nx, Ny);
fprintf('- %d Rx antennas (%d per user, %d users)\n', Nr, Nr_per_user, N_users);
fprintf('- %d OFDM subcarriers \n', N_sc);
fprintf('- Modulation order: %d\n', mod_order);
fprintf('- Cell radius: %d m, Min distance: %d m\n', cell_radius, min_distance);
fprintf('- Number of trials: %d\n\n', num_trials);
