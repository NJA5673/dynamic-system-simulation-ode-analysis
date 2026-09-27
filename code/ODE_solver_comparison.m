% ODE_Solver_Comparison.m
% Dynamic System Simulation and ODE Solver Analysis
%
% This script compares multiple numerical ODE solvers on a nonlinear
% second-order dynamic system:
%
%       y'' + sin(y) = cos(3t)
%
% The system is rewritten as a first-order state-space model and solved
% using fixed-step and adaptive-step numerical integration methods.

clear; clc; close all;

% Save plots alongside this project regardless of MATLAB's current folder.
projectRoot = fileparts(fileparts(mfilename('fullpath')));
resultsDir = fullfile(projectRoot, 'results');
if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end

% --- Figure formatting ---
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

% --- State-space model ---
% State variables:
%   x1 = y
%   x2 = dy/dt
%
% Therefore:
%   x1' = x2
%   x2' = cos(3t) - sin(x1)

systemDynamics = @(t, x) [ ...
    x(2);
    cos(3*t) - sin(x(1))
];

% --- Simulation setup ---
timeSpan = [0 50];
initialState = [0 1];

% Fixed-step solver settings
stepEuler = 0.05;
stepHeun = 0.5;
stepRK4 = 0.5;
stepButcher5 = 0.5;

% Adaptive RKF45 solver settings
tolerance = 1e-3;
minStep = 1e-4;
maxStep = 0.5;

% --- Numerical simulation ---
[timeEuler, stateEuler] = Euler_Method(systemDynamics, timeSpan, initialState, stepEuler);
[timeHeun, stateHeun] = Heun_Method(systemDynamics, timeSpan, initialState, stepHeun);
[timeRK4, stateRK4] = RungeKutta4_Method(systemDynamics, timeSpan, initialState, stepRK4);
[timeButcher5, stateButcher5] = Butcher5_Method(systemDynamics, timeSpan, initialState, stepButcher5);
[timeRKF45, stateRKF45, stepRKF45] = Adaptive_RKF45_Method( ...
    systemDynamics, timeSpan, initialState, tolerance, maxStep, minStep);

% ============================================================
% Figure 1: Position comparison
% ============================================================

figure('Name', 'Position Response Comparison', 'Position', [100 100 1200 500]);

plot(timeEuler, stateEuler(:,1), 'LineWidth', 1.5); hold on;
plot(timeRK4, stateRK4(:,1), 'LineWidth', 1.5);
plot(timeHeun, stateHeun(:,1), 'LineWidth', 1.5);
plot(timeButcher5, stateButcher5(:,1), 'LineWidth', 1.5);
plot(timeRKF45, stateRKF45(:,1), 'LineWidth', 1.5);

grid on;
xlabel('$t$ (s)', 'FontSize', 14);
ylabel('$y(t)$', 'FontSize', 14);
title('Position Response of Nonlinear Dynamic System', 'FontSize', 14);
legend({'Euler', 'RK4', 'Heun', 'Butcher 5th Order', 'Adaptive RKF45'}, ...
    'Location', 'best');
set(gca, 'FontSize', 12);

% Save figure for portfolio/report use
exportgraphics(gcf, fullfile(resultsDir, 'position_response_comparison.png'), 'Resolution', 300);

% ============================================================
% Figure 2: Velocity comparison
% ============================================================

figure('Name', 'Velocity Response Comparison', 'Position', [150 150 1200 500]);

plot(timeEuler, stateEuler(:,2), 'LineWidth', 1.5); hold on;
plot(timeRK4, stateRK4(:,2), 'LineWidth', 1.5);
plot(timeHeun, stateHeun(:,2), 'LineWidth', 1.5);
plot(timeButcher5, stateButcher5(:,2), 'LineWidth', 1.5);
plot(timeRKF45, stateRKF45(:,2), 'LineWidth', 1.5);

grid on;
xlabel('$t$ (s)', 'FontSize', 14);
ylabel('$\dot{y}(t)$', 'FontSize', 14);
title('Velocity Response of Nonlinear Dynamic System', 'FontSize', 14);
legend({'Euler', 'RK4', 'Heun', 'Butcher 5th Order', 'Adaptive RKF45'}, ...
    'Location', 'best');
set(gca, 'FontSize', 12);

% Save figure for portfolio/report use
exportgraphics(gcf, fullfile(resultsDir, 'velocity_response_comparison.png'), 'Resolution', 300);

% ============================================================
% Figure 3: Adaptive RKF45 step-size behavior
% ============================================================

figure('Name', 'Adaptive Step-Size Behavior', 'Position', [200 200 1200 500]);

stepMidpoints = timeRKF45(1:end-1) + diff(timeRKF45)/2;
plot(stepMidpoints, stepRKF45, 'LineWidth', 1.5);

grid on;
xlabel('$t$ (s)', 'FontSize', 14);
ylabel('$h(t)$', 'FontSize', 14);
title('Adaptive Step-Size Selection in RKF45 Method', 'FontSize', 14);
set(gca, 'FontSize', 12);

% Save figure for portfolio/report use
exportgraphics(gcf, fullfile(resultsDir, 'adaptive_step_size_behavior.png'), 'Resolution', 300);