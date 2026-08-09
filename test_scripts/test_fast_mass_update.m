clc; clear; close all;

rng(42);

%% This script checks the implementation for fast_mass_update
% The plots show that forming the mass matrix term scales as O(r^3) as 
% discussed in the paper.

%% Setup Parameters
r_values = 10:10:240; % Range of r values to test
num_r = length(r_values);

time_fast = zeros(num_r, 1);
time_naive = zeros(num_r, 1);
errors = zeros(num_r, 1);

%% Run Benchmark Loop
fprintf('Running benchmarks...\n');
for i = 1:num_r
    r = r_values(i);
    E_r = rand(r, r^2);
    x = rand(r, 1);
    
    % Measure fast method
    tic;
    ans_fast = fast_mass_update(E_r, x, r);
    time_fast(i) = toc;
    
    % Measure naive method
    tic;
    ans_naive = E_r * (kron(eye(r), x) + kron(x, eye(r)));
    time_naive(i) = toc;
    
    % Store error to ensure mathematical equivalence
    errors(i) = norm(ans_fast - ans_naive, 'fro');
end

fprintf('Maximum difference between methods across all r: %e\n', max(errors));

%% Plot Results
f = figure('Name', 'Scaling Comparison', 'Position', [100, 100, 1200, 500], 'Color', 'w');
set(f, 'DefaultAxesFontSize', 12); 
set(f, 'DefaultAxesTickLabelInterpreter', 'latex');

% --- Subplot 1: Linear Scale ---
subplot(1, 2, 1);
plot(r_values, time_fast, '-o', 'LineWidth', 2, 'DisplayName', 'Fast Method');
hold on;
plot(r_values, time_naive, '-s', 'LineWidth', 2, 'DisplayName', 'Naive Method');
xlabel('Reduced order ($r$)', 'Interpreter', 'latex');
ylabel('Time (seconds)', 'Interpreter', 'latex');
title('Linear Scale', 'Interpreter', 'latex');
legend('Location', 'northwest', 'Interpreter', 'latex');
grid on;

% --- Subplot 2: Log-Log Scale ---
subplot(1, 2, 2);
loglog(r_values, time_fast, '-o', 'LineWidth', 2, 'DisplayName', 'Fast Method');
hold on;
loglog(r_values, time_naive, '-s', 'LineWidth', 2, 'DisplayName', 'Naive Method');

ref_r3 = (r_values.^3) * (time_fast(end) / r_values(end)^3);
ref_r4 = (r_values.^4) * (time_naive(end) / r_values(end)^4);

loglog(r_values, ref_r3, '--k', 'LineWidth', 1.5, 'DisplayName', '$\mathcal{O}(r^3)$ Reference');
loglog(r_values, ref_r4, '-.k', 'LineWidth', 1.5, 'DisplayName', '$\mathcal{O}(r^4)$ Reference');

xlabel('Reduced order ($r$)', 'Interpreter', 'latex');
ylabel('Time (seconds)', 'Interpreter', 'latex');
title('Log-Log Scale (Complexity Verification)', 'Interpreter', 'latex');
legend('Location', 'northwest', 'Interpreter', 'latex');
grid on;
