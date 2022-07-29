%%


addpath( genpath( "./dataBuilders"))

%%
% precession parameters
s_ph   = 0.6;      % phase precession controller 1 --> perfect  -1 --> none
sig    = 8;        % width for place cell Gaussian

% *network* set up
nAsC = 10;       % number of representative place cells
nBaC = 40;       % number of background noisy cells % 80%

% timing set up 
T = 1000;            % Final time in seconds T=200 gives ~ 4 runs
bin_width = 0.01/14;  % in seconds // 0.001 sec = 1 ms


[spM, X, Times, t, position] = build_spike_data_noise(s_ph, ...
                                  nAsC, ...
                                  nBaC, ...
                                  T, ...
                                  bin_width, ...
                                  sig); 

t = t - t(1);
Times = Times - t(1);

%%
%pos vs time graph
figure;
plot(t,position)
xlabel('time')
ylabel('position')

%%
figure;




%spikes
for i = 1:size(X,1)
    
    plot(t,X(i,:),"|")
    
    drawnow
    pause(0.1)
    %hold on;
end
xlabel('time')
ylabel('spikes per bin')


%%
figure;

%looking at theta for 1 spikes per bin graph
plot(t,X(1,:),"|")
xlim([Times(1)+5,Times(1)+15])
ylim([0.99,1.01])
xlabel('time')
ylabel('spikes per bin')

h = xline(0:1/7:250);
set(h,'color','g','linewidth',0.05)



%%%for j = 1:size(spM,1)
    %%5for i = 1:size(X,1)
        
        %%plot(spM(j,:),position)
      
        %drawnow
        %pause(0.5)
    %end
%end
%%

num_delays = 10; %we like 10
tau = 1; %dont increase tau
[delayed_data, new_time] = build_delay_coordinates(X', ...
                                                  t, ...
                                                  num_delays, ...
                                                  tau);






%fast svd
[vecs, vals, modes] = svd(delayed_data, 'econ'); %
vals = diag( vals );
%peakin = zeros(length(peak_time_series(vecs(i),1:200:t(end),num_delays)),3);
% =========================================================================
%graphs for different eigenseries
%%
for i = 1:3
    figure
    
    
    plot(new_time, vecs(:,i), 'Color', [0 0 0 0.5])
    %[~,peakin] = peak_time_series(vecs(i),1:200:t(end),num_delays);

    xlim([Times(1)+5,Times(1)+15])
    
    xlabel('time')
    ylabel(sprintf('eigen time series %d', i))
    title(sprintf('bin width = %0.6f, num delays = %d, s_{ph} = %0.2f', bin_width, num_delays, s_ph))
    
    %adding lines for theta preession
    ll = 0 : (1/7) : T;
    line([1 1]'.* ll, ylim, 'color', [0.5 0.5 0.5 0.5], 'linewidth', 1);


end
%%
interval = (1/7)/bin_width;
[~,peakin] = peak_time_series(vecs(:,1),1:interval:length(t),num_delays);
ysum = 0;
xsum = 0;
x_theta = (1:numel(peakin))*(1/7);
closetotimes = zeros(1,length(Times));

for i = 1:numel(Times)
    %figure
    %plot(x_theta,peakin/interval);
    
    %xsum = xsum + 1:numel(peakin)*(1/7) == (Times(i)) : ((1:numel(peakin))*(1/7)  == (Times(i)+20));
    %ysum = ysum + peakin(  ((1:numel(peakin))*(1/7)) == (Times(i)) : ((1:numel(peakin))*(1/7)) == (Times(i)+20))./interval;
    closetotimes(i)= find( abs(x_theta - Times(i)) < 1e-10);
    

    %xlim([Times(i),Times(i)+20])
    %xlabel("time")
    %ylabel("'rightness' relative to theta period")
    %title(sprintf('indices with interval %d, bin width = %0.6f, num delays = %d, s_{ph} = %0.2f', i, bin_width, num_delays, s_ph))
end


figure
sumoflines = 0;
for i = 1:numel(Times) - 1
    hold on
    %plot(peakin(closetotimes(i):closetotimes(i)+(20*7))/interval);
    sumoflines = sumoflines + peakin(closetotimes(i):closetotimes(i)+(20*7))/interval;

    
end

trackrun = (0:140)*(1/7);
avgoflines = sumoflines/(length(closetotimes)-1);
plot(trackrun, avgoflines);
xlabel('run time')
ylabel("'rightness' of spikes")
title(sprintf('average of runs, bin width = %0.6f, num delays = %d, s_{ph} = %0.2f', bin_width, num_delays, s_ph))

line([1 1]' .* [0.4 0.6]*20, ylim, 'Color', 'r', 'LineWidth', 2)

%function to find peak in each theta period


%hyper parameter tuning
    %binsize
    %length of run
    %num delays
    %size of delays


%fourier = fft(vecs)


% xdata --> x indices for range [0.4 0.6]
% ydata --> values calculated at the indices for xdata


xindices = (floor(0.4*length(trackrun)) : ceil(0.6*length(trackrun)));
xdata  = trackrun(xindices);
ydata = avgoflines(xindices);
tbl = table(xdata(:), ydata(:));

mdl = fitlm( tbl );

intercept = mdl.Coefficients.Estimate(1);
slope = mdl.Coefficients.Estimate(2);

fit_func = @(x) intercept + slope .* x;

plot( xdata, fit_func( xdata ), 'Color', 'b', 'LineWidth', 2)

pval = mdl.Coefficients.pValue(2);
R_sq = 1 - mdl.SSE / mdl.SST;

text(max(xlim), ...
    1, ...
    sprintf('slope = %0.4f', slope), ...
    'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'right', ...
    'Interpreter', 'Latex', 'FontSize', 18)
text(max(xlim), ...
    0.9, ...
    sprintf('p value = %0.2E', pval), ...
    'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'right', ...
    'Interpreter', 'Latex', 'FontSize', 18)
text(max(xlim), ...
    0.8, ...
    sprintf('$R^2$ = %0.4f', R_sq), ...
    'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'right', ...
    'Interpreter', 'Latex', 'FontSize', 18)

ylim([0 1])


%correlation
[rho, pval] = corr(position(1:end-num_delays)',vecs(:,1));

[rho2, pval2] = corr(position(1+num_delays:end)',vecs(:,1));
