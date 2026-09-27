function [t, y] = Euler_Method(f, tspan, y0, h)
% EULER_METHOD
% First-order explicit Euler solver for systems of ODEs.
%
% This method provides a simple numerical approximation for the time
% evolution of dynamic systems governed by y' = f(t, y). While computationally
% inexpensive, it is known to accumulate error rapidly and is primarily used
% for baseline comparison with higher-order methods.
%
% Inputs:
%   f      - function handle defining the system dynamics (column output)
%   tspan  - time interval [t0 tf]
%   y0     - initial condition vector
%   h      - time step
%
% Outputs:
%   t      - time vector
%   y      - numerical solution matrix

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

    % Compute derivative
    derivative = f(currentTime, currentState);

    % Forward Euler update
    nextState = currentState + h * derivative;

    y(i + 1, :) = nextState.';
end

end