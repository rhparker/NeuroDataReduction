function [m, indices] = peak_time_series(x,t_indices,num_delays)
    %only positiveee
    x = abs(x);
    %shifting the indicessss
    new_indices = t_indices - num_delays;
    %initializing
    m = zeros(1,length(new_indices)-1);
    indices = zeros(1,length(new_indices)-1);

    %reshifting
    new_indices(1) = 1;
    
    

    %assigning
    for i = 1:length(new_indices)-1
        [m(i),indices(i)] = max(x(new_indices(i):new_indices(i+1)),[],'all');

    end



end