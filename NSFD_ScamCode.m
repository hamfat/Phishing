% MATLAB Code for Finite Difference Solution of Scam Dynamics System
clc; clear all; close all;

% Define beta values for variation
beta_values = [0.005, 0.01, 0.02, 0.03];

% Define other fixed parameters
%beta = 0.0125; %0.02;    % Rate at which susceptible individuals become victims
gamma = 0.060339; %0.01;    % Rate at which victims recover
sigma = 0.023856; %0.05;   % Rate at which recovered individuals become susceptible again
delta = 0.033928;  %0.02;   % Growth rate of active scammers
mu = 0.002901; %0.01;       % Scam/Natural exit rate of active scammers
lambda = 0.037933;  %0.1;  % Rate at which active scammers are removed due to law enforcement
psi = 0.000003; %0.005; %Rate at which some victims joined active scammers

% Time parameters
tspan = [0 1500];
dt = 0.1;
t = tspan(1):dt:tspan(2); % Using a fixed dt for consistent time vectors across simulations
N = length(t);

% Pre-allocate cell arrays to store results for each beta
num_simulations = length(beta_values);
S_results = cell(1, num_simulations);
V_results = cell(1, num_simulations);
R_results = cell(1, num_simulations);
As_results = cell(1, num_simulations);
Rs_results = cell(1, num_simulations);

% --- Main Simulation Loop ---
for k = 1:num_simulations
    % Set the beta for the current simulation
    beta = beta_values(k);
    fprintf('Running simulation for beta = %.4f\n', beta);

    % The denominator function
    dt1 = (1 - exp(-1*(sigma + mu + gamma + lambda + psi)*dt))/(sigma + mu + gamma + lambda + psi);
    
    % Initial conditions (normalized to population percentages)
    S0 = 20000;      % Initial susceptible population
    V0 = 500;       % Initial victim population
    R0 = 0.0;       % Initial recovered population
    As0 = 200;       % Initial active scammers
    Rs0 = 0.0;       % Initial removed scammers

    % Initialize solution arrays for the current simulation
    S = zeros(1, N);
    V = zeros(1, N);
    R = zeros(1, N);
    As = zeros(1, N);
    Rs = zeros(1, N);

    % Set initial conditions
    S(1) = S0;
    V(1) = V0;
    R(1) = R0;
    As(1) = As0;
    Rs(1) = Rs0;

    % Finite difference solution for the current beta
    for i = 1:(N-1)
        % Update using NSFD method
        S(i+1) = (S(i) + dt1*sigma*R(i)) / (1 + dt1*beta*As(i));
        V(i+1) = (V(i) + dt1*beta*S(i)*As(i)) / (1 + dt1*gamma + dt1*psi);
        R(i+1) = (R(i) + dt1*gamma*V(i)) / (1 + dt1*sigma);
        As(i+1) = ((As(i) + dt1*delta*As(i)) + dt1*psi*V(i)) / (1 + dt1*(mu + lambda));
        Rs(i+1) = Rs(i) + dt1*lambda*As(i);
    end
    
    % Store the results of the current simulation
    S_results{k} = S;
    V_results{k} = V;
    R_results{k} = R;
    As_results{k} = As;
    Rs_results{k} = Rs;

    % Check conservation and print R0 for the current simulation
    N1 = S + V + R + As + Rs; % total population
    R0_val = (beta * psi)*mean(N1) ./ ((gamma + psi)*(mu + lambda - delta));  %reproduction number (rate)
    fprintf('R0 for beta=%.4f: %.2f\n\n', beta, R0_val);
end

% --- Plotting Section ---
%figure('Name');

% Define colors and legend labels for the plot
colors = ['b', 'r', 'k', 'm', 'g']; % Add more if needed
legend_labels = cell(1, num_simulations);
for k = 1:num_simulations
    legend_labels{k} = sprintf('\\fontsize{13}\\beta = %.3f', beta_values(k));
end

% Plot Susceptible Population (S)
figure(1)
hold on;
for k = 1:num_simulations
    plot(t, S_results{k}, 'Color', colors(k), 'LineWidth', 1.5);
end
hold off;
title('Susceptible Population (S)');
xlabel('Time (days)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Population', 'Interpreter', 'latex', 'FontSize', 14);
legend(legend_labels);
grid on;
box on;


% Plot Victim Population (V)
figure(2)
hold on;
for k = 1:num_simulations
    plot(t, V_results{k}, 'Color', colors(k), 'LineWidth', 1.5);
end
hold off;
title('Victim Population (V)');
xlabel('Time (days)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Population', 'Interpreter', 'latex', 'FontSize', 14);
legend(legend_labels);
grid on;
box on;


% Plot Recovered Population (R)
figure(3)
hold on;
for k = 1:num_simulations
    plot(t, R_results{k}, 'Color', colors(k), 'LineWidth', 1.5);
end
hold off;
title('Recovered Population (R)');
xlabel('Time (days)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Population', 'Interpreter', 'latex', 'FontSize', 14);
legend(legend_labels);
grid on;
box on;


% Plot Active Scammers (As)
figure(4)
hold on;
for k = 1:num_simulations
    plot(t, As_results{k}, 'Color', colors(k), 'LineWidth', 1.5);
end
hold off;
title('Active Scammers (As)');
xlabel('Time (days)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Population', 'Interpreter', 'latex', 'FontSize', 14);
legend(legend_labels);
grid on;
box on;


% Plot Removed Scammers (Rs)
figure(5)
hold on;
for k = 1:num_simulations
    plot(t, Rs_results{k}, 'Color', colors(k), 'LineWidth', 1.5);
end
hold off;
title('Removed Scammers (Rs)');
xlabel('Time (days)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Population', 'Interpreter', 'latex', 'FontSize', 14);
legend(legend_labels);
grid on;
box on;
