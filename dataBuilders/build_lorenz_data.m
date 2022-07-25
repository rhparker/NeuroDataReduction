function [X,keep_time] = build_lorenz_data(noise_multiplier, t_span, dt)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% suggestion for starting/ending time
t0 = t_span(1);
tf = t_span(2);

% suggestion for temporal step (interpolation step not integration)

% for removing transients
burn_in = 10;

% parameters for chaos
params.sigma = 10;
params.rho = 28;
params.beta = 8/3;

% parameter for noise
params.eta = noise_multiplier;

% default initial position
X0 = [1 1 1];

% full time vector including transients
full_time = t0 : dt : (tf + burn_in);
% get ode solution from matlab
% sol = ode45(@(t,X) lorenz_rhs(X,params), [full_time(1) full_time(end)], X0);
X = noisy_lorenz( full_time, X0, params );


% interpolate at time points
% X = deval( sol, full_time);

% find time points past burn in
keep_time_idxs = full_time >= (t0+burn_in);

% keep only those times and set first time to t0
keep_time = full_time( keep_time_idxs ) - burn_in;

% keep state variable at same time points
X = X( :, keep_time_idxs )';

%

function sde_sol = noisy_lorenz(full_time,X0,params)
    dt = full_time(2) - full_time(1);
    
    noise_mat = params.eta * randn(numel(X0), numel(full_time));
    
    sde_sol = nan( numel(X0), numel(full_time) );
    sde_sol(:,1) = X0;
    for iter = 2:numel( full_time )
        dX = lorenz_rhs( sde_sol(:,iter-1), params );
        
        sde_sol(:,iter) =sde_sol(:,iter-1) + ...
                            dX * dt + ...
                            noise_mat(:,iter) * sqrt(dt);
    end
    
end


function dX = lorenz_rhs(X, params)
    % rhs of lorenz model from:
    %       https://en.wikipedia.org/wiki/Lorenz_system

    % extract state variables
    x = X(1);
    y = X(2);
    z = X(3);

    % extract parameters
    rho = params.rho;
    sigma = params.sigma;
    beta = params.beta;

    % governing eqs
    dx = sigma * (y - x);
    dy = x .* (rho - z) - y;
    dz = x .* y - beta * z;

    % return state changes
    dX = [dx; dy; dz];
end

fprintf('\n\n** Finished Lorenz simulation**\n\n')
end

