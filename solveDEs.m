clear
clc

params = readmatrix('lhs_params.csv');
n = size(params,1);
results = zeros(n, 2);

for i = 1:n
    tic
    p = params(i,:);
    beta = p(1); gamma = p(2); sigma = p(3);
    delta = p(4); mu = p(5); lambda = p(6); psi = p(7);

    y0 = [20000; 500; 0; 200; 0];
    tspan = [0 1520];

    ode = @(t, y) [
        -beta*y(1)*y(4) + sigma*y(3);
         beta*y(1)*y(4) - gamma*y(2) - psi*y(2);
         gamma*y(2) - sigma*y(3);
         delta*y(4) - mu*y(4) - lambda*y(4) + psi*y(2);
         lambda*y(4)
    ];

    try
        [t, y] = ode15s(ode, tspan, y0);
        cum_As = trapz(t, y(:,4));
        cum_V = trapz(t, y(:,2));
    catch
        cum_As = NaN; cum_V = NaN;
    end

    results(i,:) = [cum_As, cum_V];

    toc
end

writematrix([params, results], 'lhs_results.csv');
