addpath( genpath('dataBuilders') );

% generate spike data
bin_size = 0.005/14;
assembly_size = 5;

[spikeTimes,Y,t,T,position] = build_spike_data_noise(0.9, assembly_size, 0, 5000, bin_size, 8);
T = T - T(1);

%% delay coordinates

% find starting times for runs
t_inds = zeros(size(t));
for ind = 1:length(t)
    [~,t_inds(ind)] = min(abs(T - t(ind)));
end

% first run
run_length = 20;
Y1 = Y(:, t_inds(1):t_inds(1)+run_length/bin_size);
T1 = T(t_inds(1):t_inds(1)+run_length/bin_size);

% average of runs
Yavg = 0*Y1;
for ind = 1:length(t)-1
    Yavg = Yavg + Y(:, t_inds(ind):t_inds(ind)+run_length/bin_size);
end
Yavg = Yavg / (length(t)-1);

% % all runs
% [delayData, newTimes] = build_delay_coordinates(Y',T,1,2);

% % first run
% [delayData, newTimes] = build_delay_coordinates(Y1',T1,1,2);

% average of all runs
[delayData, newTimes] = build_delay_coordinates(Yavg',T1,1,1);

smoothDelayData = smoothdata( delayData, 'gaussian', 200);

%% SVD
% [U,S,V] = svd(delayData,'econ');
[U,S,V] = svd(smoothDelayData,'econ');

%% plot interesting stuff

% theta rhythm
thetas = 1:(1/7):200;
thetay = ones(size(thetas));
thetay(1:2:end) = 0;

figure;

% plot eigen-timeseries
eig_ind = 1;            % index of eigen-timeseries
subplot(2,1,1);
plot(newTimes, abs( U(:,eig_ind) ) );
hold on;
% superimpose theta rhythm
stairs( thetas, max(abs( U(:,eig_ind) ) )*thetay, 'r' );
xlim([t(1) t(1)+20]);

% plot spike data
subplot(2,1,2);
hold on;
for ind = 1:assembly_size
    plot( newTimes, ind+(delayData(:,ind)>0 ), 'k');
end
% superimpose theta rhythm
stairs( thetas, 6*thetay, 'r' );
xlim([t(1) t(1)+20]);

% plot animal position
% subplot(3,1,3);
% plot(T, position );
% xlim([t(1) t(1)+20]);

%%

[pks,locs] = findpeaks(abs( U(:,eig_ind) ), 'MinPeakHeight',0.005);
pktimes = newTimes(locs);
figure;
plot(mod(pktimes, 1/7))



