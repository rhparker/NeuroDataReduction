function [delayed_data, new_time] = build_delay_coordinates(x, ...
                                                  time, ...
                                                  num_delays, ...
                                                  step)
% =========================================================================
% On input: x --> data matrix size: [ num observ  x  observation dimension]
%                 i.e., each column is a timeseries observation
%           time --> time at sampled data points
%           num_delays --> amount of history to append 
%           step --> number of samples between delays (typically 1)
%
% On exit:  delayed_data --> matrix
%           let row(j) be the jth row of matrix x:
%
%           |row( n )  row( n  - step)  ... row( n  - step*num_delays)|
%           |row(n+1)  row(n+1 - step)  ... row(n+1 - step*num_delays)|
%           |...                                                      |
%           |row(N)    row( N  - step)  ... row( N  - step*num_delays)|
%           new_time --> sampled time after truncation from initial delays
% =========================================================================

decay_rate = 0;

delayed_data = x;

% add time delays one at a time
for ind = 1:num_delays
    delayed_data = [ delayed_data circshift(x,ind*step) ];    
end

% chop off top
delayed_data = delayed_data(step*num_delays+1:end,:);
new_time = time(step*num_delays+1:end,:);

end