% Nonlinear_Dynamics_and_Chaos_Simulation.m
% Analysis of Nonlinear Oscillations and Chaos Using the Duffing Equation
%
% This script investigates nonlinear dynamic behavior in a forced Duffing
% oscillator. The system exhibits periodic, multi-periodic, and chaotic
% responses depending on forcing amplitude.
%
% The analysis includes:
% - Time-domain response
% - Phase-space trajectories
% - Poincare maps
% - Strange attractor visualization
% - Flowing Poincare map comparison

clear; clc; close all;

% Save plots alongside this project regardless of MATLAB's current folder.
projectRoot = fileparts(fileparts(mfilename('fullpath')));
resultsDir = fullfile(projectRoot, 'results');
if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end

%% Figure formatting

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

%% Duffing oscillator parameters

delta = 0.3;         % damping coefficient
alpha = -1;          % linear stiffness coefficient
beta  = 1;           % nonlinear stiffness coefficient
omega = 1.2;         % forcing frequency

forcingPeriod = 2*pi / omega;

stepsPerPeriod = 260; % Practical default; increase for convergence study
timeStep = forcingPeriod / stepsPerPeriod;
poincareStep = stepsPerPeriod;

initialState = [1; 0];

%% Case definitions

cases = struct( ...
    'label', { ...
        'Period-1 Oscillation', ...
        'Weak Period-2 Oscillation', ...
        'Strong Period-2 Oscillation', ...
        'Period-4 Oscillation', ...
        'Period-5 Oscillation'}, ...
    'gamma', {0.20, 0.28, 0.65, 0.29, 0.37}, ...
    'timeEnd', {200, 200, 200, 200, 200}, ...
    'fileTag', { ...
        'period1', ...
        'weak_period2', ...
        'strong_period2', ...
        'period4', ...
        'period5'} ...
);

%% Run periodic and multi-periodic cases

for caseIndex = 1:numel(cases)

    forcingAmplitude = cases(caseIndex).gamma;
    timeEnd = cases(caseIndex).timeEnd;

    systemDynamics = @(t, y) Duffing_RHS(t, y, delta, alpha, beta, forcingAmplitude, omega);

    [time, state] = RungeKutta4_Method(systemDynamics, [0 timeEnd], initialState.', timeStep);

    position = state(:, 1);
    velocity = state(:, 2);

    transientStartIndex = find(time >= 40, 1, 'first');

    phasePosition = position(transientStartIndex:end);
    phaseVelocity = velocity(transientStartIndex:end);

    sampleTimes = ceil(40/forcingPeriod)*forcingPeriod:forcingPeriod:timeEnd;
    poincarePosition = interp1(time, position, sampleTimes);
    poincareVelocity = interp1(time, velocity, sampleTimes);

    [xLimits, vLimits] = ComputePhaseLimits(phasePosition, phaseVelocity);

    figure('Name', cases(caseIndex).label, 'Position', [100 100 1300 450]);
    tiledlayout(1, 3);

    nexttile;
    plot(time, position, 'LineWidth', 1.1);
    grid on;
    xlim([0 timeEnd]);
    xlabel('Time, $t$');
    ylabel('$x(t)$');
    title(sprintf('Response, $\\gamma = %.2f$', forcingAmplitude));

    nexttile;
    plot(phasePosition, phaseVelocity, 'LineWidth', 1.0);
    grid on;
    xlim(xLimits);
    ylim(vLimits);
    xlabel('$x$');
    ylabel('$\dot{x}$');
    title('Phase Space');

    nexttile;
    plot(poincarePosition, poincareVelocity, '.', 'MarkerSize', 7);
    grid on;
    xlim(xLimits);
    ylim(vLimits);
    xlabel('$x$');
    ylabel('$\dot{x}$');
    title('Poincare Map');

    exportgraphics(gcf, fullfile(resultsDir, sprintf('duffing_%s_response_phase_poincare.png', ...
        cases(caseIndex).fileTag)), 'Resolution', 300);

end

%% Chaotic case

chaoticGamma = 0.50;
longTimeEnd = 4000; % Practical default; increase after checking memory

chaoticDynamics = @(t, y) Duffing_RHS(t, y, delta, alpha, beta, chaoticGamma, omega);

[longTime, longState] = RungeKutta4_Method(chaoticDynamics, [0 longTimeEnd], initialState.', timeStep);

longPosition = longState(:, 1);
longVelocity = longState(:, 2);

responseEndIndex = find(longTime <= 300, 1, 'last');
transientStartIndex = find(longTime >= 40, 1, 'first');

responseTime = longTime(1:responseEndIndex);
responsePosition = longPosition(1:responseEndIndex);

phasePosition = longPosition(transientStartIndex:responseEndIndex);
phaseVelocity = longVelocity(transientStartIndex:responseEndIndex);

sampleTimes = ceil(40/forcingPeriod)*forcingPeriod:forcingPeriod:longTimeEnd;
poincarePosition = interp1(longTime, longPosition, sampleTimes);
poincareVelocity = interp1(longTime, longVelocity, sampleTimes);

[xLimits, vLimits] = ComputePhaseLimits(phasePosition, phaseVelocity);

%% Figure: Chaotic response, phase space, and Poincare map

figure('Name', 'Chaotic Duffing Response', 'Position', [100 100 1300 450]);
tiledlayout(1, 3);

nexttile;
plot(responseTime, responsePosition, 'LineWidth', 1.0);
grid on;
xlim([0 300]);
xlabel('Time, $t$');
ylabel('$x(t)$');
title('Chaotic Response, $\gamma = 0.50$');

nexttile;
plot(phasePosition, phaseVelocity, 'LineWidth', 0.5);
grid on;
xlim(xLimits);
ylim(vLimits);
xlabel('$x$');
ylabel('$\dot{x}$');
title('Phase Space, $40 \leq t \leq 300$');

nexttile;
plot(poincarePosition, poincareVelocity, '.', 'MarkerSize', 3);
grid on;
xlim(xLimits);
ylim(vLimits);
xlabel('$x$');
ylabel('$\dot{x}$');
title(sprintf('Poincare Map, $40 \leq t \leq %.0f$', longTimeEnd));

exportgraphics(gcf, fullfile(resultsDir, 'duffing_chaotic_response_phase_poincare.png'), 'Resolution', 300);

%% Figure: Strange attractor only

figure('Name', 'Duffing Strange Attractor', 'Position', [100 100 700 600]);

plot(poincarePosition, poincareVelocity, '.', 'MarkerSize', 3);
grid on;
xlim(xLimits);
ylim(vLimits);
xlabel('$x$', 'FontSize', 14);
ylabel('$\dot{x}$', 'FontSize', 14);
title('Duffing Oscillator Poincare Map: Strange Attractor', 'FontSize', 14);
set(gca, 'FontSize', 12);

exportgraphics(gcf, fullfile(resultsDir, 'duffing_strange_attractor.png'), 'Resolution', 300);

%% Flowing Poincare maps

phaseOffsets = 0.05:0.05:1.0;

figure('Name', 'Flowing Poincare Maps', 'Position', [100 100 1300 900]);

tileLayout = tiledlayout(4, 5);
title(tileLayout, 'Flowing Poincare Maps for Chaotic Duffing Oscillator, $\gamma = 0.50$');

for j = 1:numel(phaseOffsets)

    phaseOffset = phaseOffsets(j);
    startTime = ceil(40/forcingPeriod)*forcingPeriod + phaseOffset*forcingPeriod;
    shiftedTimes = startTime:forcingPeriod:longTimeEnd;
    shiftedPoincarePosition = interp1(longTime, longPosition, shiftedTimes);
    shiftedPoincareVelocity = interp1(longTime, longVelocity, shiftedTimes);

    nexttile;
    plot(shiftedPoincarePosition, shiftedPoincareVelocity, '.', 'MarkerSize', 3);
    grid on;
    xlim(xLimits);
    ylim(vLimits);
    xlabel('$x$');
    ylabel('$\dot{x}$');
    title(sprintf('$n = %.2f$', phaseOffset));

end

exportgraphics(gcf, fullfile(resultsDir, 'duffing_flowing_poincare_maps.png'), 'Resolution', 300);

%% Animated Poincare map frames

figure('Name', 'Animated Poincare Map', 'Position', [100 100 700 600]);

frames(numel(phaseOffsets)) = struct('cdata', [], 'colormap', []);

for j = 1:numel(phaseOffsets)

    phaseOffset = phaseOffsets(j);
    startTime = ceil(40/forcingPeriod)*forcingPeriod + phaseOffset*forcingPeriod;
    shiftedTimes = startTime:forcingPeriod:longTimeEnd;
    shiftedPoincarePosition = interp1(longTime, longPosition, shiftedTimes);
    shiftedPoincareVelocity = interp1(longTime, longVelocity, shiftedTimes);

    clf;
    plot(shiftedPoincarePosition, shiftedPoincareVelocity, '.', 'MarkerSize', 4);
    grid on;
    xlim(xLimits);
    ylim(vLimits);
    xlabel('$x$');
    ylabel('$\dot{x}$');
    title(sprintf('Animated Poincare Map, phase offset = %.2f', phaseOffset));

    drawnow;
    frames(j) = getframe(gcf);

end

movie(gcf, frames, 1, 4);

%% Local functions

function dydt = Duffing_RHS(t, y, delta, alpha, beta, gamma, omega)
% DUFFING_RHS
% First-order state-space representation of the forced Duffing oscillator.
%
% Original equation:
%   x'' + delta*x' + alpha*x + beta*x^3 = gamma*cos(omega*t)
%
% State vector:
%   y = [x; xdot]

position = y(1);
velocity = y(2);

acceleration = -delta*velocity - alpha*position - beta*position^3 + gamma*cos(omega*t);

dydt = [velocity; acceleration];

end

function [xLimits, vLimits] = ComputePhaseLimits(position, velocity)
% COMPUTEPHASELIMITS
% Computes padded axis limits for phase-space and Poincare plots.

xMin = min(position);
xMax = max(position);

vMin = min(velocity);
vMax = max(velocity);

padding = 0.1 * max([xMax - xMin, vMax - vMin]);

xLimits = [xMin - padding, xMax + padding];
vLimits = [vMin - padding, vMax + padding];

end