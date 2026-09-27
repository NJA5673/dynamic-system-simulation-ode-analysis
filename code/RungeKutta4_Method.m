function [t, y] = RungeKutta4_Method(f, tspan, y0, h)
% RUNGEKUTTA4_METHOD
% Classical 4th-order Runge-Kutta solver for systems of ODEs.
%
% This method provides a high-accuracy numerical solution by evaluating
% the system dynamics at multiple intermediate points within each time step.
% It achieves a global error of O(h^4), making it significantly more accurate
% than lower-order methods such as Euler and Heun.
%
% Inputs:
%   f      - function handle defining system dynamics (column output)
%   tspan  - time interval [t0 tf]
%   y0     - initial condition vector
%   h      - time step
%
% Outputs:
%   t      - time vector
%   y      - numerical solution matrix (rows = time steps)

% --- Time setup ---
t0 = tspan(1);
tf = tspan(2);

numSteps = round((tf - t0) / h);
t = linspace(t0, tf, numSteps + 1).';
h = t(2) - t(1);

% --- State initialization ---
numStates = numel(y0);
y = zeros(numSteps + 1, numStates);
y(1, :) = y0;

% --- Time integration loop ---
for i = 1:numSteps

    currentTime = t(i);
    currentState = y(i, :).';

    % Evaluate slopes at intermediate points
    k1 = f(currentTime, currentState);
    k2 = f(currentTime + h/2, currentState + (h/2)*k1);
    k3 = f(currentTime + h/2, currentState + (h/2)*k2);
    k4 = f(currentTime + h,   currentState + h*k3);

    % Weighted average update
    nextState = currentState + (h/6) * (k1 + 2*k2 + 2*k3 + k4);

    y(i + 1, :) = nextState.';
end

end