% Triple_Pendulum_Dynamics_Analysis.m
% Multi-Body Dynamic System Simulation Using Numerical ODE Solvers
%
% This script models the nonlinear motion of a triple pendulum system and
% compares multiple numerical integration methods for chaotic multi-body
% dynamics. The analysis evaluates solver agreement, computational
% efficiency, and adaptive step-size behavior.

clear; clc; close all;

% Save plots alongside this project regardless of MATLAB's current folder.
projectRoot = fileparts(fileparts(mfilename('fullpath')));
resultsDir = fullfile(projectRoot, 'results');
if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end

%% Figure formatting

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

%% Physical parameters

m1 = 3;
m2 = 2;
m3 = 1;

l1 = 1;
l2 = 1;
l3 = 1;

g = 9.81;

%% Initial conditions

theta0  = pi/2;
phi0    = pi;
psi0    = pi/4;

dtheta0 = 0;
dphi0   = 1;
dpsi0   = -1;

initialState = [theta0; phi0; psi0; dtheta0; dphi0; dpsi0];

timeStart = 0;
timeEnd = 15;
timeSpan = [timeStart timeEnd];

%% System dynamics

systemDynamics = @(t, y) TriplePendulum_RHS(t, y, m1, m2, m3, l1, l2, l3, g);

%% Solver settings

stepHeun = 0.01;
stepRK4 = 0.02;
stepButcher5 = 0.04;

tolerance = 1e-5;
maxStep = 0.1;
minStep = 1e-4;

%% Run simulations

[timeHeun, stateHeun] = Heun_Method(systemDynamics, timeSpan, initialState.', stepHeun);
[timeRK4, stateRK4] = RungeKutta4_Method(systemDynamics, timeSpan, initialState.', stepRK4);
[timeButcher5, stateButcher5] = Butcher5_Method(systemDynamics, timeSpan, initialState.', stepButcher5);

[timeRKF45, stateRKF45, stepRKF45] = Adaptive_RKF45_Method( ...
    systemDynamics, timeSpan, initialState, tolerance, maxStep, minStep);

stepTimeRKF45 = timeRKF45(1:end-1) + diff(timeRKF45) / 2;

%% Ensure row-based state matrices

stateHeun = NormalizeStateMatrix(stateHeun);
stateRK4 = NormalizeStateMatrix(stateRK4);
stateButcher5 = NormalizeStateMatrix(stateButcher5);
stateRKF45 = NormalizeStateMatrix(stateRKF45);

%% Extract angular positions

thetaHeun = stateHeun(:, 1);
phiHeun   = stateHeun(:, 2);
psiHeun   = stateHeun(:, 3);

thetaRK4 = stateRK4(:, 1);
phiRK4   = stateRK4(:, 2);
psiRK4   = stateRK4(:, 3);

thetaButcher5 = stateButcher5(:, 1);
phiButcher5   = stateButcher5(:, 2);
psiButcher5   = stateButcher5(:, 3);

thetaRKF45 = stateRKF45(:, 1);
phiRKF45   = stateRKF45(:, 2);
psiRKF45   = stateRKF45(:, 3);

%% Extract angular velocities

dthetaHeun = stateHeun(:, 4);
dphiHeun   = stateHeun(:, 5);
dpsiHeun   = stateHeun(:, 6);

dthetaRK4 = stateRK4(:, 4);
dphiRK4   = stateRK4(:, 5);
dpsiRK4   = stateRK4(:, 6);

dthetaButcher5 = stateButcher5(:, 4);
dphiButcher5   = stateButcher5(:, 5);
dpsiButcher5   = stateButcher5(:, 6);

dthetaRKF45 = stateRKF45(:, 4);
dphiRKF45   = stateRKF45(:, 5);
dpsiRKF45   = stateRKF45(:, 6);

%% Figure 1: Angular displacement response

figure('Name', 'Triple Pendulum Angle Response', 'Position', [100 100 1200 500]);

plot(timeRK4, thetaRK4, 'LineWidth', 1.1); hold on;
plot(timeHeun, thetaHeun, '--', 'LineWidth', 1.1);
plot(timeButcher5, thetaButcher5, ':', 'LineWidth', 1.1);
plot(timeRKF45, thetaRKF45, '-.', 'LineWidth', 1.1);

plot(timeRK4, phiRK4, 'LineWidth', 1.1);
plot(timeHeun, phiHeun, '--', 'LineWidth', 1.1);
plot(timeButcher5, phiButcher5, ':', 'LineWidth', 1.1);
plot(timeRKF45, phiRKF45, '-.', 'LineWidth', 1.1);

plot(timeRK4, psiRK4, 'LineWidth', 1.1);
plot(timeHeun, psiHeun, '--', 'LineWidth', 1.1);
plot(timeButcher5, psiButcher5, ':', 'LineWidth', 1.1);
plot(timeRKF45, psiRKF45, '-.', 'LineWidth', 1.1);

xline(12, 'k:', '12 s', 'LabelVerticalAlignment', 'bottom', ...
    'LabelOrientation', 'horizontal');

grid on;
xlabel('Time, $t$ (s)', 'FontSize', 14);
ylabel('Angle (rad)', 'FontSize', 14);
title('Triple Pendulum Angular Response', 'FontSize', 14);

legend({ ...
    'RK4 $\theta$', 'Heun $\theta$', 'Butcher5 $\theta$', 'RKF45 $\theta$', ...
    'RK4 $\phi$', 'Heun $\phi$', 'Butcher5 $\phi$', 'RKF45 $\phi$', ...
    'RK4 $\psi$', 'Heun $\psi$', 'Butcher5 $\psi$', 'RKF45 $\psi$'}, ...
    'Location', 'bestoutside');

set(gca, 'FontSize', 12);

exportgraphics(gcf, fullfile(resultsDir, 'triple_pendulum_angle_response.png'), 'Resolution', 300);

%% Figure 2: Angular velocity response

figure('Name', 'Triple Pendulum Angular Velocity Response', 'Position', [100 100 1200 500]);

plot(timeRK4, dthetaRK4, 'LineWidth', 1.1); hold on;
plot(timeHeun, dthetaHeun, '--', 'LineWidth', 1.1);
plot(timeButcher5, dthetaButcher5, ':', 'LineWidth', 1.1);
plot(timeRKF45, dthetaRKF45, '-.', 'LineWidth', 1.1);

plot(timeRK4, dphiRK4, 'LineWidth', 1.1);
plot(timeHeun, dphiHeun, '--', 'LineWidth', 1.1);
plot(timeButcher5, dphiButcher5, ':', 'LineWidth', 1.1);
plot(timeRKF45, dphiRKF45, '-.', 'LineWidth', 1.1);

plot(timeRK4, dpsiRK4, 'LineWidth', 1.1);
plot(timeHeun, dpsiHeun, '--', 'LineWidth', 1.1);
plot(timeButcher5, dpsiButcher5, ':', 'LineWidth', 1.1);
plot(timeRKF45, dpsiRKF45, '-.', 'LineWidth', 1.1);

xline(12, 'k:', '12 s', 'LabelVerticalAlignment', 'bottom', ...
    'LabelOrientation', 'horizontal');

grid on;
xlabel('Time, $t$ (s)', 'FontSize', 14);
ylabel('Angular Velocity (rad/s)', 'FontSize', 14);
title('Triple Pendulum Angular Velocity Response', 'FontSize', 14);

legend({ ...
    'RK4 $\dot{\theta}$', 'Heun $\dot{\theta}$', 'Butcher5 $\dot{\theta}$', 'RKF45 $\dot{\theta}$', ...
    'RK4 $\dot{\phi}$', 'Heun $\dot{\phi}$', 'Butcher5 $\dot{\phi}$', 'RKF45 $\dot{\phi}$', ...
    'RK4 $\dot{\psi}$', 'Heun $\dot{\psi}$', 'Butcher5 $\dot{\psi}$', 'RKF45 $\dot{\psi}$'}, ...
    'Location', 'bestoutside');

set(gca, 'FontSize', 12);

exportgraphics(gcf, fullfile(resultsDir, 'triple_pendulum_angular_velocity_response.png'), 'Resolution', 300);

%% Computational timing comparison

fprintf('\nComputational timing comparison on [0, 15] seconds:\n');

tic;
Heun_Method(systemDynamics, timeSpan, initialState.', stepHeun);
cpuHeun = toc;

tic;
RungeKutta4_Method(systemDynamics, timeSpan, initialState.', stepRK4);
cpuRK4 = toc;

tic;
Butcher5_Method(systemDynamics, timeSpan, initialState.', stepButcher5);
cpuButcher5 = toc;

tic;
Adaptive_RKF45_Method(systemDynamics, timeSpan, initialState, tolerance, maxStep, minStep);
cpuRKF45 = toc;

fprintf('Heun Method:       %.4f s\n', cpuHeun);
fprintf('RK4 Method:        %.4f s\n', cpuRK4);
fprintf('Butcher5 Method:   %.4f s\n', cpuButcher5);
fprintf('Adaptive RKF45:    %.4f s\n', cpuRKF45);

%% Figure 3: Adaptive step-size behavior

figure('Name', 'Triple Pendulum Adaptive Step Size', 'Position', [100 100 1200 500]);

plot(stepTimeRKF45, stepRKF45, 'LineWidth', 1.2); hold on;
yline(stepRK4, 'r--', 'RK4 fixed step', 'LineWidth', 1.2);

grid on;
xlabel('Time, $t$ (s)', 'FontSize', 14);
ylabel('Step Size, $h$', 'FontSize', 14);
title('Adaptive RKF45 Step Size Compared with Fixed RK4 Step Size', 'FontSize', 14);
legend({'Adaptive RKF45 step size', 'RK4 fixed step size'}, 'Location', 'best');
set(gca, 'FontSize', 12);

exportgraphics(gcf, fullfile(resultsDir, 'triple_pendulum_adaptive_step_size.png'), 'Resolution', 300);

%% RKF45 step-size statistics

totalRKF45Steps = numel(stepRKF45);
smallerThanRK4 = sum(stepRKF45 < stepRK4);
percentSmaller = 100 * smallerThanRK4 / totalRKF45Steps;

fprintf('\nAdaptive RKF45 step-size summary:\n');
fprintf('Total RKF45 steps: %d\n', totalRKF45Steps);
fprintf('Steps smaller than RK4 step size %.4f: %.2f%%\n', stepRK4, percentSmaller);

%% Local functions

function dydt = TriplePendulum_RHS(~, y, m1, m2, m3, l1, l2, l3, g)
% TRIPLEPENDULUM_RHS
% Defines the first-order nonlinear state-space model for the triple
% pendulum system.
%
% State vector:
%   y = [theta; phi; psi; dtheta; dphi; dpsi]

theta  = y(1);
phi    = y(2);
psi    = y(3);

dtheta = y(4);
dphi   = y(5);
dpsi   = y(6);

massMatrix = [ ...
    (m1 + m2 + m3)*l1,                  (m2 + m3)*l2*cos(theta - phi),  m3*l3*cos(theta - psi);
    (m2 + m3)*l1*cos(theta - phi),      (m2 + m3)*l2,                   m3*l3*cos(phi - psi);
    l1*cos(theta - psi),                l2*cos(phi - psi),              l3];

forcingVector = [ ...
    -((m2 + m3)*l2*sin(theta - phi)*dphi^2 + ...
      m3*l3*sin(theta - psi)*dpsi^2 + ...
      (m1 + m2 + m3)*g*sin(theta));

     (m2 + m3)*l1*sin(theta - phi)*dtheta^2 - ...
      m3*l3*sin(phi - psi)*dpsi^2 - ...
      (m2 + m3)*g*sin(phi);

      l1*sin(theta - psi)*dtheta^2 + ...
      l2*sin(phi - psi)*dphi^2 - ...
      g*sin(psi)];

angularAcceleration = massMatrix \ forcingVector;

dydt = [dtheta; dphi; dpsi; angularAcceleration];

end

function stateMatrix = NormalizeStateMatrix(stateMatrix)
% NORMALIZESTATEMATRIX
% Ensures numerical solver output is stored as rows = time steps and
% columns = state variables.

stateMatrix = squeeze(stateMatrix);

if size(stateMatrix, 2) ~= 6 && size(stateMatrix, 1) == 6
    stateMatrix = stateMatrix.';
end

end