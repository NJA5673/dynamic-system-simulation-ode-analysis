function [t, y, stepHistory] = Adaptive_RKF45_Method(f, tspan, y0, tol, hmax, hmin)
% ADAPTIVE_RKF45_METHOD
% Adaptive Runge-Kutta-Fehlberg (RKF45) solver for systems of ODEs.
%
% This method dynamically adjusts the step size based on local truncation
% error estimates by comparing embedded 4th- and 5th-order solutions.
% It provides an efficient balance between accuracy and computational cost.
%
% Inputs:
%   f      - function handle defining system dynamics (column output)
%   tspan  - time interval [t0 tf]
%   y0     - initial condition vector
%   tol    - desired error tolerance
%   hmax   - maximum allowable step size
%   hmin   - minimum allowable step size
%
% Outputs:
%   t            - time vector (non-uniform spacing)
%   y            - numerical solution matrix
%   stepHistory  - adaptive step sizes used during integration

% --- Initialization ---
t0 = tspan(1);
tf = tspan(2);

t = t0;
y = y0(:).';
stepHistory = [];

h = hmax;
i = 1;
maxAttempts = 1000000;
attempts = 0;

% --- Adaptive integration loop ---
while t(i) < tf
    attempts = attempts + 1;
    if attempts > maxAttempts
        error('Adaptive_RKF45:TooManySteps', 'Integration did not finish.');
    end

    % Ensure final step lands exactly on tf
    hi = min(h, tf - t(i));
    ti = t(i);
    yi = y(i, :).';

    % --- RKF45 slope evaluations ---
    k1 = hi * f(ti, yi);
    k2 = hi * f(ti + hi/4, yi + k1/4);
    k3 = hi * f(ti + 3*hi/8, yi + 3*k1/32 + 9*k2/32);
    k4 = hi * f(ti + 12*hi/13, yi + 1932*k1/2197 - 7200*k2/2197 + 7296*k3/2197);
    k5 = hi * f(ti + hi, yi + 439*k1/216 - 8*k2 + 3680*k3/513 - 845*k4/4104);
    k6 = hi * f(ti + hi/2, yi - 8*k1/27 + 2*k2 - 3544*k3/2565 + ...
                                 1859*k4/4104 - 11*k5/40);

    % --- Embedded solutions ---
    y4 = yi + 25*k1/216 + 1408*k3/2565 + 2197*k4/4104 - k5/5;
    y5 = yi + 16*k1/135 + 6656*k3/12825 + 28561*k4/56430 ...
              - 9*k5/50 + 2*k6/55;

    % --- Error estimate ---
    errorEstimate = max(abs(y5 - y4));

    % --- Step size adjustment factor ---
    if errorEstimate == 0
        scaleFactor = 4.0;
    else
        scaleFactor = 0.84 * (tol / errorEstimate)^(1/5);
    end

    % --- Accept step if within tolerance ---
    if errorEstimate <= tol
        t(i+1, 1) = ti + hi;
        y(i+1, :) = y5.';
        stepHistory(end+1, 1) = hi;
        i = i + 1;
    end

    if errorEstimate > tol && hi <= hmin * (1 + 1e-12)
        error('Adaptive_RKF45:MinimumStep', ...
            'Tolerance could not be met at the minimum step size.');
    end

    % --- Update step size ---
    if scaleFactor <= 0.1
        h = 0.1 * hi;
    elseif scaleFactor >= 4.0
        h = 4.0 * hi;
    else
        h = scaleFactor * hi;
    end

    % Enforce bounds
    h = min(max(h, hmin), hmax);


end

end
