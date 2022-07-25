function [spM, Y] = build_spike_data(s, sig)
% =========================================================================
% Pake Melland July 8, 2022
%   Adapted for spike data from:
%       precession.m
%       Katie Hedrick, 3/29/20
%       Model of theta phase precession
% =========================================================================


%% Parameters
v = 5;       % animal velocity (cm/sec)
omega = 7;   % LFP (Hz)
% s = 0.99;     % s=-1 for no theta, s approaches 1 for perfect phase control

rpk = 10;    % peak firing rate (Hz)
% sig = 8;     % std of Gaussian place field (cm)
rad = sig*sqrt(2*log(2*rpk));  % distance at which r/rpk = .01 (fld radius, cm)

psi = @(d) pi*(1-d);  % phase(displacement), d=(x(t)-c_i)/radius, c_i=center

ell = 100;     % track length (cm)

cvec = 40:60;  % place field center for each cell
Ncells=length(cvec);

dt = 1e-4;  % timestep (sec)


%% Set the theta modulation function
sgn = @(x) double(x>=0);

% Normalize using average firing rate across field
x = linspace(-rad,rad,1e5); t = x/v; d = x/rad; dx = x(2)-x(1);
int = dx * sum(sgn(-s+cos(2*pi*omega*t - psi(d))) .* exp(-x.^2/(2*sig^2)));
R = rpk/(2*rad) * int;
Rbar = rpk*sig*sqrt(2*pi) / (2*rad);
cth = Rbar/R;

thmod = @(t,d) cth * sgn(-s + cos(2*pi*omega*t - psi(d)));

%% Plot LFP, theta modulation, and firing rate of one cell
c0 = cvec(ceil(Ncells/2));

dx = v*dt;  % spatial step (cm)
xx = 0:dx:ell; tt = xx/v;
dd = (xx-c0)/rad;  % proportional displacement from field center

LFP = cos(2*pi*omega*tt);  % local field potential
Gfld = rpk*exp(-(xx-c0).^2/(2*sig^2));  % firing rate w/o theta modulation (Hz)
r = thmod(tt,dd) .* Gfld;  % firing rate (Hz)

figure
subplot(2,1,1)
plot(xx,LFP)
ylabel('LFP')

hold on
xpks = v*(0 : 1/omega : ell/v);  % peaks of theta rhythm
for i=1:length(xpks)
    plot(xpks(i)*[1 1],[-1 1],'--','color',.7*[1 1 1]); 
end
hold off
xlim([c0-rad, c0+rad])

subplot(2,1,2)
plot(xx,r,'b',xx,Gfld,'--b')

hold on
for i=1:length(xpks)
    plot(xpks(i)*[1 1],[0 max(r)],'--','color',.7*[1 1 1])
end
hold off
ylabel('r(t) (Hz)'); xlabel('x(t) (cm)')
xlim([c0-rad, c0+rad])

%% Plot spikes of all cells
figure
for i=1:length(xpks)
    plot(xpks(i)*[1 1],[0 Ncells+1],'color',.7*[1 1 1])
    hold on
end

S = zeros(Ncells,length(xx));
for n=1:Ncells
    c=cvec(n);
    Gfld = rpk*exp(-(xx-c).^2/(2*sig^2));  % firing rate w/o theta modulation (Hz)
    r = thmod(tt, (xx-c)/rad) .* Gfld;  % firing rate (Hz) (1 x Nsteps)
    lam = r*dt;  % 1 x Nsteps, avg # spks in each timestep
    
    spk = poissrnd(lam,1,length(xx));  % # spks in each bin
    S(n,:)=spk;
    
    xs = xx(find(spk)); ee=0*xs + n;    
    line([xs;xs],[ee-.4; ee+.4],'color','k')
end

hold off
cmin=min(cvec); cmax=max(cvec); xlim([cmin-20 cmax+20])
ylabel('Cell'); xlabel('Animal Position (cm): x = vt')
title('Poisson Spikes')
ylim([.5 Ncells+.5])

P = get(gcf,'position'); P(3)=P(3)*2; set(gcf,'position',P)



%% Compute avg firing rate across each field (check eqns)
T = 2*rad/v;    % time spent within the field
R = sum(S,2)/T; % avg firing rate for each cell
figure
plot([0 Ncells],Rbar*[1 1],'k',[0 Ncells],mean(R)*[1 1],'--r',...
    1:Ncells,R,'.b')
axis tight
xlabel('Cell')
ylabel('Firing Rate (Hz)')
title('Avg Firing Rate within Place Field')
legend('Intended','Actual','Individual')


%% Compute phase of each spike
figure
clear ph_spk x_spk
ph_spk{Ncells}=[]; d_spk{Ncells}=[];

% Plot mean phase
dd = (xx-ell/2)/rad;
plot(dd,psi(dd),'k')
hold on
plot(dd,psi(dd)+pi,'--k'); plot(dd,psi(dd)-pi,'--k')
        
for n=1:Ncells
    ind = find(S(n,:));
    d_spk{n} = (xx(ind) - cvec(n))/rad;
    ph_spk{n} = mod(omega*tt(ind),1)*2*pi;
    
    % +/- 2*pi to phase to better see and quantify the trend
    dev = psi(d_spk{n}) - ph_spk{n};  % deviation from mean phase
    ind = find(dev>pi); ph_spk{n}(ind) = ph_spk{n}(ind) + 2*pi;
    ind = find(dev<-pi); ph_spk{n}(ind) = ph_spk{n}(ind) - 2*pi;
    
    plot(d_spk{n},ph_spk{n},'.','markersize',12)
end

hold off
xlabel('Proportional Displacement from Field Center')
ylabel('Phase (rad)')
xlim([-1.5 1.5])

set(gca,'ytick',-2*pi : pi/2 : 3*pi);
set(gca,'yticklabel',{'-2\pi','-3\pi/2','-\pi','-\pi/2','0',...
    '\pi/2','\pi','3\pi/2','2\pi','5\pi/2','3\pi'})
grid on

%% PM addition
% =========================================================================
% First build spM matrix as in A. Barreiro files (where are they again?)
max_num_spikes = max( sum( S,2 ) );
spM = nan(Ncells, max_num_spikes);
for j = 1 : Ncells
    ttj = tt( S(j,:) > 0);
    spM(j,1:numel(ttj)) = ttj;
end
%%

final_time = max( max( spM ) );
start_time = min( min( spM ) );
bin_width = 0.001; % in seconds
bin_edges = (start_time-bin_width) : bin_width : final_time;


if bin_edges(end) ~= final_time
    bin_edges = [bin_edges bin_edges(end)+bin_width];
end
bin_midpts = ( bin_edges(1:end-1) + bin_edges(2:end) ) ./ 2;


Y = nan( Ncells, numel( bin_midpts ) );

for row = 1 : Ncells
    Y( row , : ) = histcounts( spM( row, : ), bin_edges ) ./ 1;
end

% figure
% imagesc( Y )
% colormap( gray )
% c = colorbar;
% c.Title.String = 'firing rate (Hz) ';
% 
% figure
% plot_levels( bin_midpts, Y(1,:), false )


% Gaussian window: make longer
nl = 20;
tau  = bin_width .* (-4*nl:4*nl);
sig_wg    =  bin_width*nl/20;
win_gauss = exp( -tau.^2/ (2 * sig_wg^2) ) / ( sqrt(2*pi) * sig_wg );
win_gauss = win_gauss ./ sum( win_gauss );


fr_gauss = conv(Y(1,:),win_gauss,'same') * bin_width;
hold on
plot( bin_midpts, fr_gauss )

%% 
X = nan( size( Y ) ) ;
for row = 1 : Ncells
    X(row,:) = conv(Y(row,:), win_gauss, 'same');
end
end

