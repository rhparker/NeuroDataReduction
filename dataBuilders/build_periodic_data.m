function [signal, t] = build_periodic_data(freq_fast, ...
                                      freq_slow, ...
                                      num_cells, ...
                                      noise_multiplier, ...
                                      t_span, ...
                                      dt, ...
                                      rng_seed)
if freq_slow == 0
    warning('''freq_slow''=0. No multiplicative modulation.')
end
rng(rng_seed)

t0 = t_span(1);
tf = t_span(2);

t = (t0 : dt : tf)'; % make time a column vector

freq_fast = freq_fast(:)';
coef_fast = rand( num_cells, numel( freq_fast) ) + 1;

coef_slow = rand( 1 , num_cells) + 5;

signal = nan( numel(t), num_cells );
beta = 100;


for cell = 1 : num_cells
    % get linear combinatin of fast frequencies
    signal_curr = sum( coef_fast(cell,:) .* cos( 2 * pi * t .* freq_fast ), 2 ); 
    
    % multiply by slow rhythm
    if freq_slow ~= 0
        signal_curr = signal_curr .* atan(beta* sin( 2 * pi * t * freq_slow ) .* coef_slow( cell ));
    end
    signal(:, cell) = signal_curr;
end

signal = signal + noise_multiplier .* randn( size(signal) );

end

