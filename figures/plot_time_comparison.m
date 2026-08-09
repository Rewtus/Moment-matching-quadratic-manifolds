clc; clear; close all;

set(0, 'defaultaxesfontsize',18,'defaultaxeslinewidth',1.0,...
    'defaultlinelinewidth',2.0,'defaultpatchlinewidth',1.0,...
    'defaulttextfontsize',18,'DefaultLineMarkerSize',14);


script_dir = fileparts(mfilename('fullpath')); 
root_dir = fullfile(script_dir, '..');         
data_dir = fullfile(root_dir, 'data');
save_path = fullfile(data_dir, 'r_results.mat');

if ~exist(save_path, 'file')
    error('Data file not found at:\n%s\nPlease run the main simulation sweep script first.', save_path);
end

fprintf('Loading simulation data from: %s\n', save_path);
load(save_path);


%% Generate Computational Runtime Plots
f1=figure('Name', 'ROM Runtime Analysis: Online vs. Offline (QM vs OSR)', ...
       'Color', 'w', 'Position', [100, 100, 1200, 500]);

% --- Subplot 1: Online Execution Time ---
subplot(1, 2, 1);
semilogy(r_values, time_qm_on,'-o'); hold on;
semilogy(r_values, time_osr_on,'-.^');
grid on; grid minor;
ylim([1e-1 10]);
xlim([min(r_values)-2, max(r_values)+2]);
xlabel('Manifold Dimension (r)',  'FontWeight', 'bold','Interpreter', 'latex');
ylabel('Execution Time [s]',  'FontWeight', 'bold','Interpreter', 'latex');
title('Online Execution Time', 'FontWeight', 'bold');
legend('QM Online', 'OSR Online', 'Location', 'southeast');

% --- Subplot 2: Offline Execution Time ---
subplot(1, 2, 2);
semilogy(r_values, time_qm_off, '-o'); hold on;
semilogy(r_values, time_osr_off, '-.^');
grid on; grid minor;

% ylim([1e-1 100]);
xlim([min(r_values)-2, max(r_values)+2]);
xlabel('Manifold Dimension (r)', 'FontWeight', 'bold','Interpreter', 'latex');
ylabel('Execution Time [s]',  'FontWeight', 'bold','Interpreter', 'latex');
title('Offline Execution Time', 'FontWeight', 'bold');
legend('QM Offline', 'OSR Offline', 'Location', 'southeast');


plot_dir = fullfile(root_dir, 'paper_plots');
if ~exist(plot_dir, 'dir')
    mkdir(plot_dir);
    fprintf('Created directory: %s\n', plot_dir);
end

print(f1, fullfile(plot_dir, 'comparison_r.eps'), '-depsc');


%% 6. Print Summary Table to Command Window
fprintf('\n==================================================================================\n');
fprintf('                     RUNTIME BENCHMARK SUMMARY (QM vs OSR)\n');
fprintf('==================================================================================\n');
fprintf('%-6s | %-14s %-14s | %-14s %-14s\n', ...
    'r', 'QM Off-Time [s]', 'QM On-Time [s]', 'OSR Off-Time [s]', 'OSR On-Time [s]');
fprintf('----------------------------------------------------------------------------------\n');
for k = 1:length(r_values)
    fprintf('%-6d | %-14.4f %-14.4f | %-14.4f %-14.4f\n', ...
        r_values(k), time_qm_off(k), time_qm_on(k), time_osr_off(k), time_osr_on(k));
end
fprintf('==================================================================================\n');