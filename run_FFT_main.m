%% For computing FFT of eigentimeseries

% uncomment to run Ross's script.
% run reduction1.m


% or you can run your own script.
% This example uses the variables as Ross's script names them.  In
% particular [U, Sigma, V] = svd( data , 'econ')
%            bin_size = width for binning spikes


eig_idx = 1;                % eigen timeseries index

Fs = 1 ./ bin_size;         % samping frequency determined by binning
signal = U(:, eig_idx);     % signal to take FFT of; column of U

[P1, freq_space] = FFT_forTimeSeries(signal, Fs);
% good idead to make title -- change variables to what you name them.
% title( sprintf('num_delays = %d, bin_width = %0.3f, s_{ph} = %0.2f', ...
%                 num_delays, bin_size, s_ph) );