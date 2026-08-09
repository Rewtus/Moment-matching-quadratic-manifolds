% Loads simulation data and generates QM ROM figures for the paper.
clc; clear; close all;
rehash;
set(0, 'defaultaxesfontsize',18,'defaultaxeslinewidth',1.0,...
    'defaultlinelinewidth',2.0,'defaultpatchlinewidth',1.0,...
    'defaulttextfontsize',18,'DefaultLineMarkerSize',14);
script_dir = fileparts(mfilename('fullpath')); 
root_dir = fullfile(script_dir, '..');         
data_path = fullfile(root_dir, 'data', 'wave_eqn_results.mat');
if ~isfile(data_path)
    error('Simulation data not found. Please run scripts/main_wave.m first.');
end
load(data_path);
fprintf('Successfully loaded simulation data from: %s\n', data_path);

Y_fom = y_t';      % FOM outputs [p x T]
Y_qm  = yr_t;      % QM ROM outputs [p x T]

% Time step size
dt = t(2) - t(1);

% --- A. Instantaneous Vector 2-Norm Error Trajectory: ||y(t) - y_r(t)||_2
e_qm_t = vecnorm(Y_fom - Y_qm, 2, 1); 

% --- B. Global Relative Output 2-Norm Error: ||Y - Y_r||_2 / ||Y||_2 ---
relerr_y_qm = norm(Y_fom - Y_qm, 2) / norm(Y_fom, 2);

% --- C. Global Relative State 2-Norm Error: ||X - X_approx||_2 / ||X||_2 ---
relerr_x_qm = norm(X_fom - X_approx, 2) / norm(X_fom, 2);

% Physical spatial grid setup
n_nodes = n / 2;
x_grid = linspace(1/n_nodes, L, n_nodes)';

%% =========================================================================
%% --- FIGURE 1: Normalized HSVs ---
%% =========================================================================
f = figure();
set(f, 'DefaultAxesTickLabelInterpreter', 'latex');
semilogy(hsv_decay / hsv_decay(1), 'o'); 
title('Normalized HSVs of the FOM', 'Interpreter', 'latex'); 
ylabel('Magnitude of HSV', 'Interpreter', 'latex'); 
xlabel('Index', 'Interpreter', 'latex'); 
xlim([1, n]); 
grid on;

%% =========================================================================
%% TABLE 2. RUNTIME STATISTICS 
%% =========================================================================
fprintf('\n==================================================================\n');
fprintf('%-15s | %-16s | %-16s | %-16s\n', 'Method', 'Offline Time (s)', 'Online Time (s)', 'Output 2Norm Err');
fprintf('------------------------------------------------------------------\n');
fprintf('%-15s | %-16.4f | %-16.4f | %-16s\n', 'FOM (Full)', 0.0000, t_fom_online, 'N/A');
fprintf('%-15s | %-16.4f | %-16.4f | %-16.4e\n', 'QM ROM', t_qm_offline, t_qm_online, relerr_y_qm);
fprintf('==================================================================\n');


%% =========================================================================
%% --- FIGURE 3: Output trajectories and absolute error ---
%% =========================================================================
fig3 = figure();
subplot(2, 1, 1);
comp = 1; % Plot first component (SISO example)
plot(t, Y_fom(comp, :), '-', 'DisplayName', 'FOM'); hold on;
plot(t, Y_qm(comp, :),  '--', 'DisplayName', 'QM ROM');
title('Steady-state output response $\mathcal M(\omega(t))$', 'Interpreter', 'latex'); 
grid on; xlabel('Time ($t$)', 'Interpreter', 'latex'); ylabel('Output Magnitude', 'Interpreter', 'latex');
legend('Location', 'northwest', 'Interpreter', 'latex');
ylim([-2*1e-5 6*1e-5]);

subplot(2, 1, 2);
semilogy(t, e_qm_t, '-', 'LineWidth', 1.5, 'DisplayName', 'Abs error');
title('Absolute Error $\|\mathcal M(\omega(t)) - \widehat{\mathcal M}(\omega(t))\|_2$', 'Interpreter', 'latex'); 
grid on; grid minor;
ylim([1e-15, 1e-10]);
xlabel('Time ($t$)', 'Interpreter', 'latex'); 
ylabel('Absolute Error', 'Interpreter', 'latex');
legend('Location', 'northwest', 'Interpreter', 'latex');

%% =========================================================================
%% --- FIGURE 4: QM State Component Analysis (Displacement & Velocity) ---
%% =========================================================================
n_half = n / 2; 

% Isolate Displacement Components (First Half of States)
FOM_disp = real(X_fom(:, 1:n_half));
ROM_disp = real(X_approx(:, 1:n_half));
err_disp = abs(FOM_disp - ROM_disp);

% Isolate Velocity Components (Second Half of States)
FOM_vel  = real(X_fom(:, n_half+1:end));
ROM_vel  = real(X_approx(:, n_half+1:end));
err_vel  = abs(FOM_vel - ROM_vel);

fig4 = figure('Name', 'State Component Analysis', 'Position', [100, 100, 1500, 900]);

% Displacement Analysis
subplot(2, 3, 1); imagesc(x_grid, t, FOM_disp); title(sprintf('FOM Displacement ($n=%d$)', n), 'Interpreter', 'latex'); xlabel('Spatial position ($z$)', 'Interpreter', 'latex'); ylabel('Time ($t$)', 'Interpreter', 'latex'); axis xy; colorbar;
subplot(2, 3, 2); imagesc(x_grid, t, ROM_disp); title(sprintf('ROM Displacement ($r=%d$)', r), 'Interpreter', 'latex'); xlabel('Spatial position ($z$)', 'Interpreter', 'latex'); ylabel('Time ($t$)', 'Interpreter', 'latex'); axis xy; colorbar;
subplot(2, 3, 3); imagesc(x_grid, t, err_disp); title('Abs. Error (Displacement)', 'Interpreter', 'latex'); xlabel('Spatial position ($z$)', 'Interpreter', 'latex'); ylabel('Time ($t$)', 'Interpreter', 'latex'); axis xy; colorbar;

% Velocity Analysis
subplot(2, 3, 4); imagesc(x_grid, t, FOM_vel); title(sprintf('FOM Velocity ($n=%d$)', n), 'Interpreter', 'latex'); xlabel('Spatial position ($z$)', 'Interpreter', 'latex'); ylabel('Time ($t$)', 'Interpreter', 'latex'); axis xy; colorbar;
subplot(2, 3, 5); imagesc(x_grid, t, ROM_vel); title(sprintf('ROM Velocity ($r=%d$)', r), 'Interpreter', 'latex'); xlabel('Spatial position ($z$)', 'Interpreter', 'latex'); ylabel('Time ($t$)', 'Interpreter', 'latex'); axis xy; colorbar;
subplot(2, 3, 6); imagesc(x_grid, t, err_vel); title('Abs. Error (Velocity)', 'Interpreter', 'latex'); xlabel('Spatial position ($z$)', 'Interpreter', 'latex'); ylabel('Time ($t$)', 'Interpreter', 'latex'); axis xy; colorbar;

colormap turbo;

%% =========================================================================
%% 2-Norm error metrics summarized
%% =========================================================================
fprintf('\n===================================================\n');
fprintf('GLOBAL RELATIVE STATE 2-NORM ERROR (X):\n');
fprintf('  QM Method:        %.4e\n', relerr_x_qm);
fprintf('---------------------------------------------------\n');
fprintf('GLOBAL RELATIVE OUTPUT 2-NORM ERROR (Y):\n');
fprintf('  QM Method:        %.4e\n', relerr_y_qm);
fprintf('===================================================\n');

plot_dir = fullfile(root_dir, 'paper_plots');
if ~exist(plot_dir, 'dir')
    mkdir(plot_dir);
    fprintf('Created directory: %s\n', plot_dir);
end

fprintf('\nSaving Figures 1, 3, and 4 as .eps files to folder: %s ...\n', plot_dir);

% Save HSV plot
print(f, fullfile(plot_dir, 'wave_hsv.eps'), '-depsc');
% Save Figure 3
print(fig3, fullfile(plot_dir, 'wave_output_trajectories.eps'), '-depsc');
% Save Figure 4
print(fig4, fullfile(plot_dir, 'wave_qm_state_reconstruction.eps'), '-depsc', '-r300');

fprintf('Successfully saved Figures 1, 3, and 4 as vector .eps for the paper!\n');