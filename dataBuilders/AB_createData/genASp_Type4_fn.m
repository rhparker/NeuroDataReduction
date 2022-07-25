function [spM] = genASp_Type4_fn(nC,freq,T,lags)
% genAsSp_Type4

% "For the less precise assembly type IV, short windows of extra 
% spikes for each assembly neuron were organized in a specific 
% temporal pattern, with the exact occurrence of the extra spikes 
% within the defined time windows determined randomly on each repetition
% (cf. Figure 1A; e.g. [Friedrich et al., 2004; Euston et al., 2007; 
% Luczak et al., 2007; Peyrache et al., 2009; Adler et al., 2012]). 
% Specifically, time windows of 0.3 s with extra spikes at a Poisson 
% rate of 10 Hz were (without loss of generality) arranged in a 
% sequential order, with the time lag between these windows drawn 
% from a uniform distribution, [0 0.4] s. While this sequential ordering
% of time windows was fixed, within each window spikes were drawn at 
% random on each assembly repetition. "
%

%

% Create lags?
if (nargin<4)
    lags =  unifrnd(0,0.4,[nC-1 1]); 
    % wrt first cell
    lags = [0;cumsum(lags)];
else
    if (length(lags)==1)
        % Interpret this as a STRETCH parameter
        stretch = lags;
        lags =  unifrnd(0,0.4,[nC-1 1]); 
        % WRT first cell
        lags = [0;cumsum(lags)]*stretch;
    end
end


tref = 0.015;  %Refractory period
% Generate times: ISIs are exponentially generated w/ mean "1/freq"
% Freq is in units of Hz (or same units as T)
%isi  = exprnd(1/freq, 1, round(freq*T*1.1));
% OR
% Use a Gamma distribution
%    Note: 1st parameter is SHAPE distribution.
nu   = 4;
isi  = gamrnd(nu,1/nu/freq,1,round(freq*T*1.2));

% Ensure above tref
isi  = isi+tref;

% Now sum
Times = cumsum(isi);

% Remove extra
Times(Times>T)=[];

% This time, pattern needs to be re-generated for each instance of the
% assembly
lTimes      = length(Times);
% For each assembly activation, how many spikes for each cell
spCount     = poissrnd(10*0.3,nC,lTimes); spCount = max(spCount,1);

lpat        = max(max(spCount));
spM = nan(nC,lpat,lTimes);
for j1=1:lTimes
    spAs     = nan(nC,lpat);
    for j2=1:nC
        % Times uniformly distributed in [0,0.3]
        ns   = spCount(j2,j1);
        spAs(j2,1:ns)=sort(unifrnd(0,0.3,[1,ns]));
    end
    % Now include lags
    spAs = spAs+lags;
   
    % Correct for violations of the refractory period
    [I,J]=find(diff(spAs,[],2)<tref);
    for j2=1:length(I)
        % I(j2) is the row: the cellID
        % J(j2) is the offending entry of "diff". That means
        %    spAs(I(j2),J(j2)+1)- spAs(I(j2),J(j2)) was too small
        %    Fix this by shifting the J(j2)+1 and later spike times
        spAs(I(j2),(J(j2)+1):end)=spAs(I(j2),(J(j2)+1):end)+tref;
    end
    
    if (any(diff(spAs,[],2)<tref))
        warning('Warning: ref period violated!');
    end
    spM(:,:,j1)=spAs+Times(j1);
end

% Now sort, remove nan
spM = reshape(spM,[nC,lTimes*lpat]);
spM = sort(spM,2);  % Along 2nd dimension

% Will probably not be any to remove
spM(:,all(isnan(spM)))=[];

