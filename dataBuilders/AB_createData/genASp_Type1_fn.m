function [spM] = genASp_Type1_fn(nC,freq,T)
% genAsSp_Type1
%

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

% Same times for all cells
spM  = repmat(Times,[nC,1]);

end

