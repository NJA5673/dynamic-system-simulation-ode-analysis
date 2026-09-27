function [t, y] = Butcher5_Method(f, tspan, y0, h)
% BUTCHER5_METHOD
% Fifth-order explicit Runge-Kutta solver for nonlinear dynamic systems.
%
% Inputs:
%   f      - system of first-order ODEs as a function handle
%   tspan  - simulation interval [t0 tf]
%   y0     - initial condition vector
%   h      - fixed time step
%
% Outputs:
%   t      - time vector
%   y      - numerical state history

t0 = tspan(1);
tf = tspan(2);

numSteps = round((tf - t0) / h);
t = linspace(t0, tf, numSteps + 1).';
h = t(2) - t(1);

numStates = numel(y0);
y = zeros(numSteps + 1, numStates);
y(1, :) = y0;

for i = 1:numSteps
    currentTime = t(i);
    currentState = y(i, :).';

    k1 = f(currentTime, currentState);
    k2 = f(currentTime + h/4, currentState + (h/4)*k1);
    k3 = f(currentTime + h/4, currentState + (h/8)*k1 + (h/8)*k2);
    k4 = f(currentTime + h/2, currentState - (h/2)*k2 + h*k3);
    k5 = f(currentTime + 3*h/4, currentState + (3*h/16)*k1 + (9*h/16)*k4);
    k6 = f(currentTime + h, currentState - (3*h/7)*k1 + (2*h/7)*k2 + ...
        (12*h/7)*k3 - (12*h/7)*k4 + (8*h/7)*k5);

    nextState = currentState + (h/90) * ...
        (7*k1 + 32*k3 + 12*k4 + 32*k5 + 7*k6);

    y(i + 1, :) = nextState.';
end

end