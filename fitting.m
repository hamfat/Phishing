clc; clear;
close all;

% Load Cleaned Monthly Victim Count Data ===
data = readtable('monthly_victim_counts.csv');
data.month = datetime(data.month, 'InputFormat', 'yyyy-MM-dd');
t_obs = days(data.month - data.month(1));  % days from first month
victim_counts = data.victim_count;

% Initial Conditions and Parameter
init = [20000, 500, 0, 200, 0];  % [S, V, R, As, Rs]
%init = [2000, 5, 0, 10, 0];
p0 = [0.02, 0.01, 0.05, 0.02, 0.01, 0.1, 0.005];  % [beta, gamma, sigma, delta, mu, lambda, psi]


% Fit parameters using fmincon
obj_fun = @(params) model_error(params, init, t_obs, victim_counts);
options = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',5000);
best_params = fmincon(obj_fun, p0, [], [], [], [], ...
    [0, 0, 0, 0, 0, 0, 0], [0.1, 0.1, 0.2, 0.1, 0.1, 1, 0.1], [], options);


% === Output Estimated Parameters ===
param_names = {'beta', 'gamma', 'sigma', 'delta', 'mu', 'lambda', 'psi'};
disp('Estimated Parameters:');
for i = 1:length(param_names)
    fprintf('%-8s = %.6f\n', param_names{i}, best_params(i));
end

% === Simulate Using Fitted Parameters ===
[~, Y_fit] = ode23s(@(t, y) scammer_ode(t, y, best_params), t_obs, init);

% === Plot Results ===
figure;
plot(t_obs, victim_counts, 'ko', 'DisplayName', 'Observed Victims'); hold on;
plot(t_obs, Y_fit(:,2), 'r-', 'LineWidth', 2, 'DisplayName', 'Model Fit');
xlabel('Time (days)'); ylabel('Victim Count');
legend; title('Fitted Model to Monthly Victim Data');
grid on;


% === Plot All State Variables
figure;
state_names = {'Susceptible (S)', 'Victims (V)', 'Recovered (R)', ...
               'Active Scammers (As)', 'Removed Scammers (Rs)'};

for i = 1:5
    subplot(3,2,i);
    plot(t_obs, Y_fit(:,i), 'b-', 'LineWidth', 2); hold on;
    if i == 2
        % Overlay observed victim counts for V
        plot(t_obs, victim_counts, 'ro', 'MarkerSize', 6, ...
             'DisplayName', 'Observed Victims');
        legend('Model Fit', 'Observed Victims', 'Location', 'best');
    else
        legend('Model Fit', 'Location', 'best');
    end
    xlabel('Time (days)');
    ylabel(state_names{i});
    title(['Fitted ', state_names{i}]);
    grid on;
end

% Adjust layout
sgtitle('Fitted Model Results for All State Variables');



% === Varying betas ===
beta_values = [0.005, 0.01, 0.02, 0.03, 0.05];  % Five specific beta values
colors = lines(length(beta_values));
figure;

% Loop through each beta value
for i = 1:length(beta_values)
    params_varied = best_params;
    params_varied(1) = beta_values(i); % Set beta to current value
    [~, Y_varied] = ode23s(@(t, y) scammer_ode(t, y, params_varied), t_obs, init);

    % Plot Victim Count (V)
    plot(t_obs, Y_varied(:,2), '-', 'Color', colors(i,:), ...
         'LineWidth', 2, 'DisplayName', sprintf('\\beta = %.3f', beta_values(i)));
    hold on;
end

xlabel('Time (days)');
ylabel('Victim Count');
title('Effect of Different β Values on Victim Dynamics');
legend('Location', 'best');
grid on;




%% === Functions ===

function err = model_error(params, init, t_obs, y_obs)
    [~, Y] = ode15s(@(t, y) scammer_ode(t, y, params), t_obs, init);
    V_model = Y(:,2); 
    err = sum((V_model - y_obs).^2);
end

function dydt = scammer_ode(~, y, p)
    beta = p(1); gamma = p(2); sigma = p(3);
    delta = p(4); mu = p(5); lambda = p(6); psi = p(7);
    S = y(1); V = y(2); R = y(3); As = y(4); Rs = y(5);
    dS  = -beta * S * As + sigma * R;
    dV  = beta * S * As - gamma * V - psi * V;
    dR  = gamma * V - sigma * R;
    dAs = delta * As - mu * As - lambda * As + psi * V;
    dRs = lambda * As;
    dydt = [dS; dV; dR; dAs; dRs];
end