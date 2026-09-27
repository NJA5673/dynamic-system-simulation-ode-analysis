function [t, y] = Heun_Method(f, tspan, y0, h)
% HEUN_METHOD
% Second-order predictor-corrector solver for systems of ODEs.
%
% This method improves on Forward Euler by first predicting the next state
% using an Euler step, then correcting the estimate using the average slope
% between the beginning and end of the time interval.
%
% Inputs:
%   f      - function handle defining system dynamics (column output)
%   tspan  - simulation interval [t0 tf]
%   y0     - initial condition vector
%   h      - fixed time step
%
% Outputs:
%   t      - time vector
%   y      - numerical state history

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

% --- Predictor-corrector integration loop ---
for i = 1:numSteps
    currentTime = t(i);
    currentState = y(i, :).';

    % Predictor step: estimate next state using Euler method
    initialSlope = f(currentTime, currentState);
    predictedState = currentState + h * initialSlope;

    % Corrector step: average initial and predicted slopes
    finalSlope = f(currentTime + h, predictedState);
    nextState = currentState + (h / 2) * (initialSlope + finalSlope);

    y(i + 1, :) = nextState.';
end

end