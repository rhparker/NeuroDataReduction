clear; clc; % clear up for this section
addpath( genpath('plottingHelpers') );
addpath( genpath('dataBuilders') );
% precession parameters
s_ph   = 0.9;      % phase precession controller 1 --> perfect  -1 --> none
sig    = 8;        % width for place cell Gaussian

% *network* set up
nAsC = 20;       % number of representative place cells
nBaC = 10;       % number of background noisy cells

% timing set up 
T = 200;            % Final time in seconds T=200 gives ~ 4 runs
bin_width = 1;  % in seconds // 0.001 sec = 1 ms

[spM, X, Times, t, position] = build_spike_data_noise(s_ph, ...
                                  nAsC, ...
                                  nBaC, ...
                                  T, ...
                                  bin_width, ...
                                  sig); 
                              
figure
plot_levels(t, X(15,:)./bin_width, false, 'Color', 'r', 'LineWidth', 4 );
hold on
stairs(t, X(15,:)./bin_width, 'Color', 'b' );