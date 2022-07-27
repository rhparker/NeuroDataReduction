function [P1, freq_space] = FFT_forTimeSeries(signal,...
                                              Fs)
% =========================================================================
% Compute fft for eigentimeseries.  Adapted from:
%       https://www.mathworks.com/help/matlab/ref/fft.html
% On input: signal --> Thing to take fft of.  Most likely, a column of U
%                      where [U, Sigma, V] = svd( data, 'econ' )
%           Fs --> sampling frequency.  Most likely, 1 ./ bin_width where
%                  bin_width is the width of spike binning
% On exit: Makes frequency plot and draws line at 7 Hz (specifically for
%          this application where the theta cycle is set to 7 Hz.
%
%          P1 --> Approximate coefficients if it were perfect sum of sines
%          freq_space --> Frequency space specified by Fs and length of
%                         signal
%          Can plot outside of this function plot( freq_space, P1 )
% =========================================================================

L = numel( signal );    % number of elements in signal

% if L isn't even, we make it even and chop our signal by one
if mod(L,2) ~= 0
    L = L-1;
end
signal = signal(1:L);

% do the Fourier Transform thanks to MATLAB
X_ft = fft(signal);

% compute one-sided power
P2 = abs( X_ft ./ L );
P1 = P2(1:L/2 + 1);
P1(2:end-1) = 2 * P1(2:end-1);

% define frequency based on sampling rate and signal length
freq_space = Fs * (0 : (L/2) ) ./ L;

% draw a plot
figure
plot(freq_space, (P1))
xlim([0 10])
line([7 7], ylim, 'Color', [0 0 0 0.25], 'LineWidth', 2)
xlabel('Frequency (Hz)')
ylabel('Approximate amplitude at frequency')
set(gca, 'FontSize', 18)
set(gca, 'LineWidth', 1.5)
s_size = get(0, 'screensize');
set(gcf, 'Position', [s_size(3)/8, ...
                      s_size(4)/4, ...
                      s_size(3)/4*3, ...
                      s_size(4)/2]);
                  
end


