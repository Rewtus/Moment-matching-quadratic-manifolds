clc; clear; close all;
rng(42); % Set a seed 
script_dir = fileparts(mfilename('fullpath')); 
root_dir = fullfile(script_dir, '..');         
addpath(fullfile(root_dir, 'src'));

data_dir = fullfile(root_dir, 'data'); % Added data_dir definition
if ~exist(data_dir, 'dir')
    mkdir(data_dir);
end

%% 1. Parameters & System Construction
n = 1024; 
L = 1;
v = 1; 
[A, B, C] = cde(n, L, v); 
A = full(A); B = full(B); C = full(C);

r = 20; % Reduced order for QM ROM
r_lin = r + r^2;
p = size(B, 2);
m = size(C, 1);

%% 2. Signal Space & Interpolation Configurations
freqs = -logspace(1, 3, r/2); 
blocks = arrayfun(@(a) [0, a; -a, 0], freqs, 'UniformOutput', false);
S = blkdiag(blocks{:});
S_2 = kron(S, eye(r)) + kron(eye(r), S); 
L_1 = 0.1*ones(p, r); 
L_2 = 0.1*ones(p, r^2); 

%% 3. Compute Interpolating QM ROM (Offline)
fprintf('--- Computing Interpolating QM ROM (Offline) ---\n');
tic; 
V1 = sylvester(-A, S, B*L_1);
V2 = sylvester(-A, S_2, B*L_2);
W = V1; 
t_qm_offline = toc; 

E1r = W'*V1;
E2r = E1r\W'*V2;
A_r = E1r\(W'*A*V1);
H_r = E1r\(W'*A*V2);
B_r = E1r\(W'*B);    
C_r = C*V1;      
K_r = C*V2;

%% 4. Define Initial Conditions & Dynamics
w0 = 0.01 * rand(r, 1);
x0_fom = V1 * w0 + V2 * kron(w0, w0); 
xr0 = w0; 

z0_fom = [w0; x0_fom];
zr0_rom = [w0; xr0]; 
fast_kron = @(x) reshape(x * x.', [], 1);

M_rom = @(t, z) blkdiag(eye(r), eye(r) + fast_mass_update(E2r, z(r+1:end), r)); 
aug_dynamics_rom = @(t, z) [
    S * z(1:r); ...
    A_r * z(r+1:end) + H_r * fast_kron(z(r+1:end)) + ...
    B_r * (L_1 * z(1:r) + L_2 * fast_kron(z(1:r)))
];

aug_dynamics_fom = @(t, z) [
    S * z(1:r); ...
    A * z(r+1:end) + B * (L_1 * z(1:r) + L_2 * kron(z(1:r), z(1:r)))
];

T = 1000;
tspan = linspace(0, 1, T);
options_rom = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'Mass', M_rom, 'OutputFcn', @ode_timer);
options_fom = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'OutputFcn', @ode_timer);

%% 5. Time Integration (Online)
fprintf('\nSolving interpolating QM ROM (Online)...\n');
tic;
[t, Z_rom] = ode15s(aug_dynamics_rom, tspan, zr0_rom, options_rom);
t_qm_online = toc;

fprintf('\nSolving FOM (Online)...\n');
tic;
[~, Z_fom] = ode15s(aug_dynamics_fom, t, z0_fom, options_fom);
t_fom_online = toc;

W_traj = Z_fom(:, 1:r);
X_fom  = Z_fom(:, r+1:end);
Xr_rom = Z_rom(:, r+1:end);

%Use first 1000 snapshots for training the greedy QM model
X_train=X_fom(1:T,:); 

%% 6. Compute Greedy QM ROM (Offline)
fprintf('\n--- Running Greedy QM Selection (Offline) ---\n');
tic;
[V_gQM, V2_gQM, selected_modes] = GreedyQM(X_train', r, 1e-8, n/4); 

W_g = V_gQM;
E1_g = W_g'*V_gQM;
E2_g = E1_g\W_g'*V2_gQM;
A_g = E1_g\(W_g'*A*V_gQM);
H_g = E1_g\(W_g'*A*V2_gQM);
B_g = E1_g\(W_g'*B);    
C_g = C*V_gQM;      
K_g = C*V2_gQM;
t_greedy_offline = toc;

S_hat = V_gQM'*X_train';
X_greedy_approx = V_gQM*S_hat + V2_gQM*compressed_quadkron(S_hat);
X_greedy_lin = V_gQM*S_hat;

%% 7. Simulate Greedy ROM (Online)
zr0_g = [w0; V_gQM'*x0_fom]; 

aug_dynamics_greedy = @(t, z) [
    S * z(1:r); ...
    A_g * z(r+1:end) + H_g * compressed_quadkron(z(r+1:end)) + ...
    B_g * (L_1 * z(1:r) + L_2 * fast_kron(z(1:r)))
];

options_g = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'OutputFcn', @ode_timer);
fprintf('\nSolving greedy ROM (Online)...\n');
tic;
[t2, Z_g] = ode15s(aug_dynamics_greedy, t, zr0_g, options_g);
t_greedy_online = toc;
Xr_g = Z_g(:, r+1:end);

%% 8. Post-Processing & State/Output Reconstructions
fprintf('\nReconstructing trajectories and outputs...\n');
X_approx   = zeros(length(t), n);
X_approx_g = zeros(length(t2), n);
y_t        = (C * X_fom')';
yr_t       = zeros(m, length(t));
yg_t       = zeros(m, length(t2));

for i = 1:length(t)
    xi_r = Xr_rom(i, :)';
    X_approx(i, :) = (V1 * xi_r + V2 * kron(xi_r, xi_r))';
    yr_t(:, i) = C_r * xi_r + K_r * kron(xi_r, xi_r);
end

for i = 1:length(t2)
    xi_g = Xr_g(i, :)';
    X_approx_g(i, :) = (V_gQM * xi_g + V2_gQM * compressed_quadkron(xi_g))';
    yg_t(:, i) = C_g * xi_g + K_g * compressed_quadkron(xi_g);
end

hsv_decay=compute_hsv(A,B,C);

save_path = fullfile(data_dir, 'qm_greedy_results.mat');
save(save_path);
fprintf('\nSimulation complete. Workspace saved successfully to "%s".\n', save_path);

%% ==================================================================
%%                       HELPER FUNCTIONS
%% ==================================================================
function [V, W, I] = GreedyQM(S, r, gamma, m)
    [Phi, Sigma, Psi] = svd(S, 'econ');
    k = size(Phi, 2); 
    I = []; I_check = 1:k; V = [];           
    A_full = Sigma * Psi'; 
    
    % figure('Name', 'Singular Value Decay');
    % semilogy(diag(Sigma), 'LineWidth', 1.5); grid on;
    % title('Snapshot Matrix Singular Values');
    
    trace_AA = sum(Sigma(:).^2);
    row_norms_sq = sum(A_full.^2, 2); 
    
    for i = 1:r
        best_val = Inf; best_j = -1;
        num_candidates = min(m, length(I_check));
        for idx = 1:num_candidates
            j_cand = I_check(idx);
            sub_indices = I_check(I_check ~= j_cand);
            hat_S_cand = A_full([I, j_cand], :); 
            H_prime = compressed_quadkron(hat_S_cand);
            H_mat = H_prime * H_prime' + gamma * eye(size(H_prime, 1));
            G = A_full * H_prime'; 
            M_G = H_mat \ (G');   
            val = (trace_AA - row_norms_sq(j_cand)) - sum(sum(G(sub_indices, :) .* M_G(:, sub_indices)'));
            if val < best_val
                best_val = val; best_j = j_cand;
            end
        end
        I = [I, best_j];
        I_check(I_check == best_j) = []; 
        V = [V, Phi(:, best_j)];
    end
    W = compute_W(V, S, gamma);
end

function y = compressed_quadkron(X)
    [n, m] = size(X);
    p = n * (n + 1) / 2;
    y = zeros(p, m); 
    row_idx = 1;
    for i = 1:n
        y(row_idx : row_idx + i - 1, :) = X(1:i, :) .* X(i, :);
        row_idx = row_idx + i;
    end
end

function W = compute_W(V, S, gamma)
    H = compressed_quadkron(V'*S);
    p = size(H, 1);
    E = (S - V * V'*S); 
    W = (E * H') / (H * H' + gamma * eye(p));
end

function status = ode_timer(t, z, flag)
    persistent last_time
    status = 0;
    switch flag
        case 'init'
            last_time = tic;
        case ''
            last_time = tic;
        case 'done'
            fprintf('Integration Step Complete.\n');
    end
end

function M_update = fast_mass_update(E2r, x, r)
    M_update = zeros(r, r);
    for i = 1:r
        E_i = E2r(:, (i-1)*r + 1 : i*r); 
        M_update = M_update + x(i) * E_i; 
        M_update(:, i) = M_update(:, i) + E_i * x; 
    end
end
