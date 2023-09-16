%Prompt: write code in matlab for stochastic runge-kutta 4th order
%Certainly! Below is a MATLAB code implementation of the Stochastic Runge-Kutta 
% 4th Order (SRK4) method to solve a simple stochastic differential equation (SDE):


function X = stochastic_rk4_solver(f, g, tspan, X0, dt)
    % f: Deterministic part of the SDE
    % g: Stochastic part of the SDE
    % tspan: Time span [t_start, t_end]
    % X0: Initial condition
    % dt: Time step size
    
    % Initialize variables
    t = tspan(1):dt:tspan(2);
    N = length(t);
    dW = sqrt(dt) * randn(1, N-1);  % Generate random noise
    
    X = zeros(size(X0, 1), N);
    X(:, 1) = X0;
    
    % SRK4 integration
    for i = 1:N-1
        % Compute deterministic increments (K1, K2, K3, K4)
        K1 = f(X(:, i), t(i));
        K2 = f(X(:, i) + 0.5*dt*K1, t(i) + 0.5*dt);
        K3 = f(X(:, i) + 0.5*dt*K2, t(i) + 0.5*dt);
        K4 = f(X(:, i) + dt*K3, t(i) + dt);
        
        % Compute stochastic increments (L1, L2, L3, L4)
        L1 = g(X(:, i), t(i)) * dW(i);
        L2 = g(X(:, i) + 0.5*dt*K1, t(i) + 0.5*dt) * dW(i);
        L3 = g(X(:, i) + 0.5*dt*K2, t(i) + 0.5*dt) * dW(i);
        L4 = g(X(:, i) + dt*K3, t(i) + dt) * dW(i);
        
        % Update the state using SRK4 formula
        X(:, i+1) = X(:, i) + (dt/6) * (K1 + 2*K2 + 2*K3 + K4) + (1/sqrt(3)) * (L1 + L2 + L3 + L4);
    end
end

% Define your SDE functions f and g, initial conditions, and time span
% Example usage:
% tspan = [0, 1];
% X0 = [0];
% dt = 0.01;
% f = @(X, t) ...  % Define your deterministic part
% g = @(X, t) ...  % Define your stochastic part
% X = stochastic_rk4_solver(f, g, tspan, X0, dt);