clc; clear; close all;

%% 1. Environment & Path Setup
script_dir = fileparts(mfilename('fullpath')); 
root_dir = fullfile(script_dir, '..');         
addpath(fullfile(root_dir, 'src'));
data_dir = fullfile(root_dir, 'data');
if ~exist(data_dir, 'dir'), mkdir(data_dir); end

%% 2. Define Parametric Sweep Values
r_values = 4:2:26; %Values chosen such that (r+r^2<n)
num_runs = length(r_values);

% Preallocate runtime metrics (Offline vs Online)
time_qm_off  = zeros(num_runs, 1);
time_qm_on   = zeros(num_runs, 1);
time_osr_off = zeros(num_runs, 1);
time_osr_on  = zeros(num_runs, 1);

%% 3. Generate Fixed FOM Parameters
n = 1024; L = 1; v = 1; 
[A, B, C] = cde(n, L, v);
A = full(A); B = full(B); C = full(C);
p = size(B, 2); 
m = size(C, 1);

% Define common simulation configurations
tspan = linspace(0, 1, 1000);
options_osr = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
fast_kron = @(x) reshape(x * x.', [], 1);

%% 4. Execute the Sweep Loop
for k = 1:num_runs
    r = r_values(k);
    r_lin = r + r^2;
    fprintf('\n==================================================\n');
    fprintf('  Evaluating ROM Order: r = %d (r_lin = %d)\n', r, r_lin);
    fprintf('==================================================\n');
    
    % Set seed inside loop
    rng(42); 
    
    % Signal space configuration
    freqs = -logspace(1, 3, r/2); 
    s = [freqs(:); -freqs(:)];
    sum_matrix = s + s';
    pairwise_sums = sum_matrix(:);
    freqs_osr_lin = sort([s; pairwise_sums]);
    freqs_osr_lin = -1j*freqs_osr_lin(1:r_lin/2);
    blocks = arrayfun(@(a) [0, a; -a, 0], freqs, 'UniformOutput', false);
    S = blkdiag(blocks{:});
    L_1 = ones(p, r); 
    L_2 = ones(p, r^2);
    w0 = 0.1 * rand(r, 1);
    
    %% --- 1. Quadratic Manifold (QM) Method ---
    fprintf(' -> Running QM...\n');
    tic;
    [V1, V2, W] = compute_qm_rom(A, B, S, L_1, L_2, r);
    E1r = W'*V1;  E2r = E1r\W'*V2;
    A_r = E1r\(W'*A*V1); H_r = E1r\(W'*A*V2); B_r = E1r\(W'*B);
    time_qm_off(k) = toc;
    
    % Online
    zr0_rom = [w0; w0]; 
    x0_fom = V1*w0+V2*kron(w0,w0); %Center manifold of the FOM
    M_rom = @(t, z) blkdiag(eye(r), eye(r) + fast_mass_update(E2r, z(r+1:end), r)); 
    aug_dynamics_rom = @(t, z) [S * z(1:r); A_r * z(r+1:end) + H_r * fast_kron(z(r+1:end)) + B_r * (L_1 * z(1:r) + L_2 * fast_kron(z(1:r)))];
    options_rom = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'Mass', M_rom);
    
    tic;
    [t, Z_rom] = ode15s(aug_dynamics_rom, tspan, zr0_rom, options_rom);
    time_qm_on(k) = toc;
    

    %% --- 2. One-Sided Rational Interpolation (OSR) ---
    fprintf(' -> Running One-Sided Rational Interpolation...\n');
    tic;
    [A_osr_r, B_osr_r, C_osr_r, V_osr] = compute_one_sided_rom(A, B, C, freqs_osr_lin, r_lin);
    time_osr_off(k) = toc;
    
    % Online
    xr0_osr = V_osr' * x0_fom;
    zr0_osr = [w0; xr0_osr];
    aug_dynamics_osr = @(t, z) [S * z(1:r); A_osr_r * z(r+1:end) + B_osr_r * (L_1 * z(1:r) + L_2 * fast_kron(z(1:r)))];
    

    tic;
    [~, Z_osr] = ode15s(aug_dynamics_osr, t, zr0_osr, options_osr);
    time_osr_on(k) = toc;
    
end

% Save updated sweep results to data_dir
save_path = fullfile(data_dir, 'r_results.mat');
save(save_path);
fprintf('Simulation complete. Data saved to %s\n', save_path);