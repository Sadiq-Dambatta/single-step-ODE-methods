% main.m
% Comparative error and convergence analysis of explicit
% single-step numerical methods for first-order ODEs
% Author: Sadiq Haruna Shehu

clear; clc; close all;

% Test problem: y' = -y, y(0) = 1, exact solution y = exp(-x)
f = @(x,y) -y;
exact = @(x) exp(-x);
x0 = 0; y0 = 1; xf = 2;

% Methods to compare
methods = {'Euler', 'Heun', 'Midpoint', 'RK4'};
h_values = [0.5, 0.25, 0.125, 0.0625, 0.03125];

fprintf('%-10s %-12s %-15s %-15s\n', 'Method', 'h', 'Error at xf', 'Order');
fprintf('-----------------------------------------------------\n');

figure; hold on; grid on;
colors = {'r-o', 'b-s', 'g-^', 'k-d'};

for m = 1:length(methods)
    method = methods{m};
    errors = zeros(size(h_values));
    for i = 1:length(h_values)
        h = h_values(i);
        N = round((xf - x0)/h);
        x = linspace(x0, xf, N+1);
        y = zeros(1, N+1);
        y(1) = y0;
        for n = 1:N
            xn = x(n); yn = y(n);
            switch method
                case 'Euler'
                    y(n+1) = yn + h*f(xn, yn);
                case 'Heun'
                    k1 = f(xn, yn);
                    k2 = f(xn+h, yn+h*k1);
                    y(n+1) = yn + h/2*(k1+k2);
                case 'Midpoint'
                    k1 = f(xn, yn);
                    k2 = f(xn+h/2, yn+h/2*k1);
                    y(n+1) = yn + h*k2;
                case 'RK4'
                    k1 = f(xn, yn);
                    k2 = f(xn+h/2, yn+h/2*k1);
                    k3 = f(xn+h/2, yn+h/2*k2);
                    k4 = f(xn+h, yn+h*k3);
                    y(n+1) = yn + h/6*(k1+2*k2+2*k3+k4);
            end
        end
        errors(i) = abs(y(end) - exact(xf));
    end
    % Convergence order
    orders = -diff(log(errors))./diff(log(h_values));
    loglog(h_values, errors, colors{m}, 'LineWidth', 1.5, 'DisplayName', method);
    for i = 1:length(h_values)
        if i == 1
            fprintf('%-10s %-12.5f %-15.6e %-15s\n', method, h_values(i), errors(i), '-');
        else
            fprintf('%-10s %-12.5f %-15.6e %-15.4f\n', method, h_values(i), errors(i), orders(i-1));
        end
    end
end

xlabel('Step size h'); ylabel('Absolute error at x = 2');
title('Convergence of Single-Step Methods for y'' = -y');
legend('show'); hold off;
