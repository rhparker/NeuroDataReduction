function [spM, Times, v, ell] = genASp_LinTrackTheta_fn(nC,freq,T,pm)
% genASp_LinTrackTheta_fn
%  Generate spikes corresponding to place fields  
%  running on a linear track 
%  Based on code from K. Hedrick: precession.m (2021)
%
%  INPUTS:
%    nC     :   # cells
%    freq   :   frequency of track run (Hz)
%    T      :   total length of spike train
%    pm     : optional parameters
%
%  PARAMETERS (units)/[range]/[default]
%    pm.v       : rat velocity (cm/sec)/[>0]/[5]
%    pm.omega   : Theta Freq (Hz)/[>0]/[7]
%    pm.s       : degree of theta modulation ()/[ [-1,1) ]/[0]
%    pm.rpk     : peak firing rate (Hz)/[>0]/[5]
%    pm.sig     : Std.dev. of Gaussian place field (cm)/[>0]/[4]
%    pm.ell     : Track length (cm)/[>0]/[100]
%    pm.tref    : Refractory period (s)/[>=0]/[0.005]
%    pm.cvec    : Place field centers (cm)/ [size=(1xnC)]/[linspace(40,60,nC)]


%% Parameters
% Default values
v = 5;       % animal velocity (cm/sec)
omega = 7;   % LFP (Hz)
s = 0;     % s=-1 for no theta, s approaches 1 for perfect phase control

rpk = 5;    % peak firing rate (Hz)
sig = 4;     % std of Gaussian place field (cm)
ell = 100;   % track length (cm)

tref = 0.005;  % Refractory period?
cvec = linspace(40,60,nC);  % place field center for each cell


%% Read any parameters passed in
if (isfield(pm,'v') && pm.v >0); v = pm.v; end
if (isfield(pm,'omega') && pm.omega >0); omega = pm.omega; end
if (isfield(pm,'s') && pm.s >= -1 && pm.s<1); s = pm.s; end
if (isfield(pm,'rpk') && pm.rpk >0); rpk = pm.rpk; end
if (isfield(pm,'sig') && pm.sig >0); sig = pm.sig; end
if (isfield(pm,'ell') && pm.ell >0); ell = pm.ell; end
if (isfield(pm,'tref') && pm.tref >=0); tref = pm.tref; end
if (isfield(pm,'cvec') && length(pm.cvec)==nC); cvec = pm.cvec; end 

% Derived
rad = sig*sqrt(2*log(2*rpk));  % distance at which r/rpk = .01 (fld radius, cm)
psi = @(d) pi*(1-d);  % phase(displacement), d=(x(t)-c_i)/radius, c_i=center

% Q: why have this be less than a ms?
dt = 1e-3;  % timestep (sec)

%% Set the theta modulation function
sgn = @(x) double(x>=0);

% Normalize using average firing rate across field
x = linspace(-rad,rad,1e5); t = x/v; d = x/rad; dx = x(2)-x(1);
int = dx * sum(sgn(-s+cos(2*pi*omega*t - psi(d))) .* exp(-x.^2/(2*sig^2)));
R = rpk/(2*rad) * int;
Rbar = rpk*sig*sqrt(2*pi) / (2*rad);
cth = Rbar/R;

thmod = @(t,d) cth * sgn(-s + cos(2*pi*omega*t - psi(d)));

% Structures we need to create spikes
dx = v*dt;  % spatial step (cm)
xx = 0:dx:ell; tt = xx/v;

% To create phase plot, if desired
c0 = cvec(ceil(nC/2));
dd = (xx-c0)/rad;  % proportional displacement from field center
xpks = v*(0 : 1/omega : ell/v);  % peaks of theta rhythm


%% Generate "run times" i.s. times when the rat will run the track

% Use a Gamma distribution
nu   = 4;
isi  = gamrnd(nu,1/nu/freq,1,round(freq*T*1.2));

% Ensure above minRunTime: rat can't be in two places at once!
minRunTime = ell/v;
isi  = isi+minRunTime;

% Now sum
Times = cumsum(isi);

% Ensure that shift times are an integer multiple of theta 
%   period
Times = round(Times*omega)/omega;

% Remove extra
Times(Times>T)=[];

% This time, pattern needs to be re-generated for each instance of the
% assembly
lTimes      = length(Times);

phFig = figure;
clear ph_spk x_spk
ph_spk{nC}=[]; d_spk{nC}=[];
% Plot mean phase
dd = (xx-ell/2)/rad;
plot(dd,psi(dd),'k')
hold on
plot(dd,psi(dd)+pi,'--k'); plot(dd,psi(dd)-pi,'--k');

spM = [];

lTimes;
Times;
for j1=1:lTimes
%% Create spikes (in terms of space)
    S = zeros(nC,length(xx));
    for n=1:nC
        c=cvec(n);
        Gfld = rpk*exp(-(xx-c).^2/(2*sig^2));  % firing rate w/o theta modulation (Hz)
        %Use theta modulation from SHIFTED time
        r = thmod(tt+Times(j1), (xx-c)/rad) .* Gfld;  % firing rate (Hz) (1 x Nsteps)
        lam = r*dt;  % 1 x Nsteps, avg # spks in each timestep

        spk = poissrnd(lam,1,length(xx));  % # spks in each bin
        S(n,:)=spk;

        %xs = xx(find(spk)); %ee=0*xs + n;    
        % line([xs;xs],[ee-.4; ee+.4],'color','k')
    end
    %% Create matrix of spike times which can be used by 
    %%     the CAD algorithm
    S    = min(S,1);   % Remove any counts >1
    nSpk = sum(S,2); maxSpk = max(nSpk);
    spMt = nan(nC,maxSpk);
    for j2=1:nC
        spMt(j2,1:nSpk(j2))= xx(find(S(j2,:)));
    end
    
    
    % Convert to units of time, add on shift
    spMt = spMt/v + Times(j1);
    
    % Remove spikes that violate a refractory period
    if (tref>0)
       [I,J]=find(diff(spMt,[],2)<tref);
       for j2=1:length(I)
           spMt(I(j2),J(j2)+1)=nan;
       end
    end
    
    %[sum(diff(spMt,[],2)<tref,2) sum(~isnan(spMt),2)]
    
    % This is to plot, phase, in case we want to check
    for n=1:nC
        % Unfortunately, we still need this
        ind = find(S(n,:));
        d_spk{n} = (xx(ind) - cvec(n))/rad;
        ph_spk{n} = mod(omega*tt(ind),1)*2*pi;

        % +/- 2*pi to phase to better see and quantify the trend
        dev = psi(d_spk{n}) - ph_spk{n};  % devition from mean phase
        ind = find(dev>pi); ph_spk{n}(ind) = ph_spk{n}(ind) + 2*pi;
        ind = find(dev<-pi); ph_spk{n}(ind) = ph_spk{n}(ind) - 2*pi;

        plot(d_spk{n},ph_spk{n},'.','markersize',12)
    end
    
    
    spM     = [spM spMt];
end

% This is for phase plot
xlabel('Proportional Displacement from Field Center')
ylabel('Phase (rad)')
xlim([-1.5 1.5])
title(sprintf('degree of theta modulation s = %0.2f', s))

% Now sort, remove nan
spM = sort(spM,2);  % Along 2nd dimension

% Remove extra nan
spM(:,all(isnan(spM)))=[];




end