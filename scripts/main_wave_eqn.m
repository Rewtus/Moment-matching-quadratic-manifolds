clc; clear; close all;
rng(42) % Set a seed
script_dir = fileparts(mfilename('fullpath')); 
root_dir = fullfile(script_dir, '..');         
addpath(fullfile(root_dir, 'src'));
data_dir = fullfile(root_dir, 'data');
if ~exist(data_dir, 'dir')
    mkdir(data_dir);
    fprintf('Created missing directory: %s\n', data_dir);
end

%% Define the FOM

n = 2048; % Number of gridpoints used in discretization
c = 1; % Speed of the wave
L = 1; % Length of the interval
gamma = 1e-8; % Low diffusion term for numerical stability

% Generate wave equation system
[A, B, C] = wave_eqn(n, c, L, gamma); 
n = 2*n; % State dimension for second-order formulation

% Dimension of the reduced systems
r = 16;

A = full(A);
B = full(B);
C = full(C);

%% Signal space parameters
% To handle MIMO case
p = size(B, 2);
m = size(C, 1);
freqs = -logspace(1, 2.5, r/2); % freqs we want to interpolate 

blocks = arrayfun(@(a) [0, a; -a, 0], freqs, 'UniformOutput', false);
S = blkdiag(blocks{:});
S_2 = (kron(S, eye(r)) + kron(eye(r), S)); 
L_1 = ones(p, r); % Simple choice in SISO case
L_2 = ones(p, r^2); % Can also pick conjugate directions 

%% =========================================================================
%% Compute ROM Matrices (OFFLINE)
%% =========================================================================
% QM Method 
fprintf('Computing QM ROM...\n');
tic; 
[V1, V2, W] = compute_qm_rom(A, B, S, L_1, L_2, r);
E1r = W'*V1;
E2r = E1r\(W'*V2);
A_r = E1r\(W'*A*V1);
H_r = E1r\(W'*A*V2);
B_r = E1r\(W'*B);    
C_r = C*V1;      
K_r = C*V2;
t_qm_offline = toc; 
fprintf('QM completed.');


% Compute HSVs
hsv_decay = compute_hsv(A, B, C);

%% Initialize the ode solver
w0 = 0.1 * rand(r, 1);
x0_fom = V1 * w0 + V2 * kron(w0, w0);

% Initial Conditions for ROMs
xr0 = w0; 

% Combined initial conditions for the ODE solvers [signal; state]
z0_fom = [w0; x0_fom];
zr0_rom = [w0; xr0];

%% =========================================================================
%% Simulate ROMs and FOM to check online performance 
%% =========================================================================

fast_kron = @(x) reshape(x * x.', [], 1);
% Optimized O(r^3) Mass matrix computation
M_rom = @(t, z) blkdiag(eye(r), eye(r) + fast_mass_update(E2r, z(r+1:end), r)); 

% Optimized QM ROM Dynamics 
aug_dynamics_rom = @(t, z) [
    S * z(1:r); ...
    A_r * z(r+1:end) + H_r * fast_kron(z(r+1:end)) + ...
    B_r * (L_1 * z(1:r) + L_2 * fast_kron(z(1:r)))
];

% FOM Dynamics (baseline)
aug_dynamics_fom = @(t, z) [
    S * z(1:r); ...
    A * z(r+1:end) + B * (L_1 * z(1:r) + L_2 * kron(z(1:r), z(1:r)))
];

% Simulation Span
tspan = linspace(0, 2, 1000);
options_rom = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'Mass', M_rom);
options_fom = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

% Compute online times 
fprintf('Solving QM ROM ...\n');
tic;
[t, Z_rom] = ode15s(aug_dynamics_rom, tspan, zr0_rom, options_rom);
t_qm_online = toc;


fprintf('Solving FOM ...\n');
tic;
[~, Z_fom] = ode15s(aug_dynamics_fom, t, z0_fom, options_fom);
t_fom_online = toc;

% Extract trajectories
W_traj      = Z_fom(:, 1:r);
X_fom       = Z_fom(:, r+1:end);
Xr_rom      = Z_rom(:, r+1:end);

% Compute the state and output reconstructions
X_approx = zeros(length(t), n);
y_t = (C * X_fom')';
yr_t = zeros(m, length(t));

for i = 1:length(t)
    xi_r = Xr_rom(i, :)';
    X_approx(i, :) = (V1 * xi_r + V2 * kron(xi_r, xi_r))';
    yr_t(:, i) = C_r * xi_r + K_r * kron(xi_r, xi_r);
end

if ~exist('data', 'dir')
    mkdir('data');
    fprintf('Created missing directory: data/\n');
end
save_path = fullfile(data_dir, 'wave_eqn_results.mat');
save(save_path);
fprintf('Simulation complete. Data saved to %s\n', save_path);