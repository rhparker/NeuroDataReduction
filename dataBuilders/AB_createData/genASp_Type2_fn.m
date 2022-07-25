function [spM] = genASp_Type2_fn(nC,freq,T,lags)
% genAsSp_Type2
%
%  "For assembly type II, spikes follow a precise sequential pattern across the set 
%  of assembly neurons on each instance of activation (Lee and Wilson, 2002; 
%   Diba and Buzsa ?ki, 2007). Time lags between spikes
%   were drawn from a uniform distribution [0 0.1] s, and then fixed for
%   each occurrence"
%

% Create lags?
if (nargin<4)
    lags =  unifrnd(0,0.1,[nC-1 1]); 
    % Order by increasing lag
    lags = [0;cumsum(lags)];
else
    if (length(lags)==1)
        % Interpret this as a STRETCH parameter
        stretch = lags;
        lags =  unifrnd(0,0.1,[nC-1 1]); 
        % Order by increasing lag
        lags = [0;cumsum(lags)]*stretch;
    end
end

tref = 0.015;  %Refractory period
% Generate "Type 1" assembly spikes: precisely times sequences

% Generate times: ISIs are exponentially generated w/ frequency "1/freq"
isi  = exprnd(freq, 1, round(freq*T*1.1));

% Ensure above tref
isi  = max(isi,tref);

% Now sum
Times = cumsum(isi);

% Remove extra
Times(Times>T)=[];

% Now add to the pattern we created earlier
spM  = repmat(Times,[nC,1])+repmat(lags,[1 length(Times)]);

end

