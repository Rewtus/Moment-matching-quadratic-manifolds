clc; clear; close all;

set(0, 'defaultaxesfontsize',18,'defaultaxeslinewidth',1.0,...
    'defaultlinelinewidth',2.0,'defaultpatchlinewidth',1.0,...
    'defaulttextfontsize',18,'DefaultLineMarkerSize',14);


script_dir = fileparts(mfilename('fullpath')); 
root_dir = fullfile(script_dir, '..');  
data_file = fullfile(root_dir, 'data', 'qm_greedy_results.mat');

if ~exist(data_file, 'file')
    error('Data file "%s" not found! Run the simulation script first.', data_file);
end
load(data_file);

%Establish Spatial Grid
xi = cos(pi * (0:n-1)' / n); 
x_grid = (L/2) * (xi + 1);
[X, T] = meshgrid(x_grid, t);

%%  Print Execution Runtime Table
fprintf('\n==================================================================\n');
fprintf('%-15s | %-20s | %-20s\n', 'Method', 'Offline Time (s)', 'Online Time (s)');
fprintf('------------------------------------------------------------------\n');
fprintf('%-15s | %-20.4f | %-20.4f\n', 'FOM (Full)', 0.0000, t_fom_online);
fprintf('%-15s | %-20.4f | %-20.4f\n', 'QM ROM', t_qm_offline, t_qm_online);
fprintf('%-15s | %-20.4f | %-20.4f\n', 'QM Greedy ROM', t_greedy_offline, t_greedy_online);
fprintf('==================================================================\n');

%%  Print Summary Accuracy Metrics
relerr = norm(X_fom - X_approx, 'fro') / norm(X_fom, 'fro');
relerr_greedy = norm(X_fom - X_approx_g, 'fro') / norm(X_fom, 'fro');

global_error_proj = norm(X_train' - X_greedy_approx, 'fro') / norm(X_fom', 'fro');
global_error_lin = norm(X_train' - X_greedy_lin, 'fro') / norm(X_fom', 'fro');

fprintf('\n===================================================\n');
if exist('selected_modes', 'var')
    fprintf('Greedy QM Selected SVD Modes: %s\n', mat2str(selected_modes));
    fprintf('---------------------------------------------------\n');
end
fprintf('--- Global Horizon [0, T] Errors ---\n');
fprintf('Global Relative Error (QM ROM):        %.4e\n', relerr);
fprintf('Global Relative Error (Greedy QM ROM): %.4e\n', relerr_greedy);
fprintf('Global Projective Error (Quad):        %.4e (%.2f%%)\n', global_error_proj, global_error_proj * 100);
fprintf('Global Projective Error (Linear):      %.4e (%.2f%%)\n', global_error_lin, global_error_lin * 100);
fprintf('===================================================\n');

%% Common Plot Styling Settings
lblX = 'Spatial Coordinates ($z$)';
lblY = 'Time ($t$)';

%% Spatiotemporal State and Error Comparison 
f1 = figure('Position', [50, 100, 1800, 700]);

% FOM State
subplot(2,3,1); 
pcolor(X, T, real(X_fom)); 
shading interp; axis xy; colorbar;
title(sprintf('FOM State ($n=%d$)', n), 'Interpreter', 'latex');
xlabel(lblX, 'Interpreter', 'latex'); ylabel(lblY, 'Interpreter', 'latex');

% QM ROM State
subplot(2,3,2); 
pcolor(X, T, real(X_approx)); 
shading interp; axis xy; colorbar;
title(sprintf('Interpolating QM ROM (r=%d)', r), 'Interpreter', 'latex'); 
xlabel(lblX, 'Interpreter', 'latex'); ylabel(lblY, 'Interpreter', 'latex');

% Greedy QM ROM State
subplot(2,3,3); 
pcolor(X, T, real(X_approx_g)); 
shading interp; axis xy; colorbar;
title(sprintf('Greedy QM ROM (r=%d)', r), 'Interpreter', 'latex'); 
xlabel(lblX, 'Interpreter', 'latex'); ylabel(lblY, 'Interpreter', 'latex');

sing_vals = svd(X_fom);

%Singular value decay of snapshot matrix
subplot(2,3,4); 
semilogy(sing_vals / sing_vals(1), 's-', 'LineWidth', 1.2, 'MarkerSize', 4); 
title('Normalized Snapshot SVD Decay', 'Interpreter', 'latex'); 
ylabel('Singular Value Magnitude', 'Interpreter', 'latex'); 
xlabel('Index', 'Interpreter', 'latex'); 
grid on;

% Abs. Error (QM ROM)
subplot(2,3,5); 
pcolor(X, T, abs(X_fom - X_approx)); 
shading interp; axis xy; colorbar;
title(sprintf('Abs. error QM ROM (r=%d)', r), 'Interpreter', 'latex'); 
xlabel(lblX, 'Interpreter', 'latex'); ylabel(lblY, 'Interpreter', 'latex');

% Abs. Error (Greedy QM ROM)
subplot(2,3,6); 
pcolor(X, T, abs(X_fom - X_approx_g)); 
shading interp; axis xy; colorbar;
title(sprintf('Abs. error gQM ROM (r=%d)', r), 'Interpreter', 'latex'); 
xlabel(lblX, 'Interpreter', 'latex'); ylabel(lblY, 'Interpreter', 'latex');

colormap turbo;


%% Output trajectories (Nonlinear moment comparison)
figure('Name', 'Outputs Comparison', 'Color', 'w');
semilogy(t, abs(y_t(:)), 'r', 'LineWidth', 1.5); hold on;
semilogy(t, abs(yr_t(:)), '--b', 'LineWidth', 1.5);

if exist('t2', 'var') && exist('yg_t', 'var')
    semilogy(t2, abs(yg_t(:)), '-.g', 'LineWidth', 1.5);
elseif exist('yg_t', 'var')
    semilogy(t, abs(yg_t(:)), '-.g', 'LineWidth', 1.5);
end
hold off;

title('System Output Trajectories ($y(t)$)', 'Interpreter', 'latex'); 
grid on; 
xlabel('Time [s]', 'Interpreter', 'latex'); 
ylabel('Magnitude (Log Scale)', 'Interpreter', 'latex');
legend('FOM Output', 'Interpolating QM ROM', 'Greedy QM ROM', 'Location', 'Best');


plot_dir = fullfile(root_dir, 'paper_plots');
if ~exist(plot_dir, 'dir')
    mkdir(plot_dir);
    fprintf('Created directory: %s\n', plot_dir);
end

fprintf('\nSaving Figures 1, as .eps files to folder: %s ...\n', plot_dir);

print(f1, fullfile(plot_dir, 'greedy_comparison_plots.eps'), '-depsc');

fprintf('Successfully saved Figure 1 as vector .eps for the paper!\n');

