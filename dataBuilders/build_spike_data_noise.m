function [newSp, Y, Times, t, position] = build_spike_data_noise(s, ...
                                             nAsC, ...
                                             nBaC, ...
                                             T, ...
                                             bin_width, ...
                                             sig)
% =========================================================================
% Pake Melland July 8, 2022
%   Adapted for spike data from:
%       precession.m
%       Katie Hedrick, 3/29/20
%       - and -
%       makeLinTrackPFspikes
%       Andrea Barreiro, 2021
%       Model of theta phase precession
% =========================================================================

addpath( genpath('AB_creatData') );

%% Things a user should control (function inputs)
% s         --> precession controllor
% nAsC      --> number of assembly cells
% nBaC      --> number of non-assembly (Background) cells
% T         --> final time
% bin_width --> bin size for binning discrete spikes
% sig       --> place cell sigma

% most other stuff, I guess can be swept under the rug to reduce my &
% student effort.

%% Output
% newSp  ---> Matrix that gives precise spike times of each cell
%             (row) x (col) <----> (num cells) x (num spikes-ish)
% Y      ---> newSp times binned according to bin_size
%             (row) x (col) <----> (num cells) x (num binned time points)
%%
% Create Poisson process as described in RD17

% nC = 10;            % Total # of cells
nC = nAsC + nBaC;            % Total # of cells
% T  = 1200;          % Total time (s) // changed to user input
lam = 0.1;          % Avg firing rate for background process
fracVar = 0.1;      % Fraction which is variable, vs. constant rate Poisson

S = makeAR(nC,T,lam);
S = fracVar*S + lam*(1-fracVar);

% Estimate size of spike matrix
maxD   = round(lam*T*1.5);
bgSpM  = nan(nC,maxD);

% Now create spikes, using a call to 
% genNHPP:   creates a spike train for a nonhomogeneous Poisson process
for j1=1:nC
    rvec = S(j1,:); tvec = 0:T;
    ETtemp  = genNHPP(@(x)eval_pw_const(x,rvec,tvec),T,1);
    
    nEv     = length(ETtemp);
    bgSpM(j1,1:nEv)=ETtemp;   
end

% Delete unneccessary columns
bgSpM(:, all(isnan(bgSpM)))=[];

% Add on refractory period as suggested by RD17
tref  = 0.015;
maxSp = size(bgSpM,2);
bgSpM = bgSpM+ones(nC,1)*[0:maxSp-1]*tref;

if (0)
% Plot
figure;
subplot(2,1,1);hold on;
for j1=1:nC
    aus = bgSpM(j1,:); aus(isnan(aus))=[];
    plot(aus,j1*ones(size(aus)),'.','color',[0.5 0.5 0.5]);
end
subplot(2,1,2);plot(0:T,S');
end

assemb_type = 10; 
% nAsC = 6;
% Make assembly
switch assemb_type
    case 1
        % TYPE I
        freq = 1;  
        spAs   = genASp_Type1_fn(nAsC,freq,T);
    case 2
        % TYPE II
        freq = 1;
        spAs   = genASp_Type2_fn(nAsC,freq,T);
    case 3
        % TYPE III
        % The larger "stretch", the more likely cells will appear in 
        %   the indicated order
        %
        % Possible issue w/ this structure: patterns are up to 200 ms long
        %      high frequency can lead to "run-in" and violations of 
        %      tref. Need to decide what to do...
        freq    = 0.5;
        stretch = 2;
        spAs   = genASp_Type3_fn(nAsC,freq,T,stretch);
    case 4
        % TYPE IV
        % The larger "stretch", the more likely cells will appear in 
        %   the indicated order
        freq    = 0.2;
        stretch = 2;
        spAs   = genASp_Type4_fn(nAsC,freq,T,stretch); 
    case 10
        % Theta modulated place field spikes
        freq = 0.04;
        pm=[]; 
%         pm.s    = 0.99;      % theta modulation parameter
        pm.s    = s;      % theta modulation parameter // changed to user input
        pm.tref = 0.002;  % refractory period
%         pm.sig  = 8;
        pm.sig  = sig;    % changed to user input
        [spAs,Times,v,ell]   = genASp_LinTrackTheta_fn(nAsC,freq,T,pm); 
end

% Am I doing refrac right?
min(diff(spAs,[],2),[],2) 


% Embed into background noise
cellIDs     = 1:nAsC;
newSp       = embedAssembly_fn(spAs,bgSpM,cellIDs);


%% PM binning
final_time = max( max( newSp ) );
start_time = 0;

% or for brevity, just take the first spike
% start_time = min( min( newSp ) );

% bin_width = 0.001; % in seconds // changed to user input
bin_edges = (start_time-bin_width) : bin_width : final_time;

if bin_edges(end) ~= final_time
    bin_edges = [bin_edges bin_edges(end)+bin_width];
end
bin_midpts = ( bin_edges(1:end-1) + bin_edges(2:end) ) ./ 2;
t = bin_midpts;
position = zeros( 1, numel(t) );
[~,Times_start_location] = min( abs((t) - Times(:)), [], 2);
Times_end = Times(:) + ell/v + 1;
[~,Times_end_location] = min( abs((t) - Times_end(:)), [], 2);


for runs = 1:numel( Times_start_location )
    pos_start = Times_start_location(runs);
    pos_end   = Times_end_location(runs);
    position_temp = linspace(0, ell, pos_end - pos_start + 1);
    
    position( pos_start : pos_end ) = position_temp;
end

Y = nan( nC, numel( bin_midpts ) );

for row = 1 : nC
    Y( row , : ) = histcounts( newSp( row, : ), bin_edges ) ./ 1;
end

%% Save data
% Some possible naming conventions:
%
% If looking at theta procession in place fields
%   fname = sprintf('test_data_s=%4.2f_lam=%4.2f.mat',pm.s,lam);
%
% If doing a "type X" assembly
%   fname = sprintf('test_data_type=%d.mat,assemb_type);
% fname = 'my_spike_trains.mat';
% save(fname,'newSp','spAs','cellIDs');

figure;
%Plot all spikes
subplot(211);hold on;
for j1=1:nC
    aus = newSp(j1,:); aus(isnan(aus))=[];
    plot(aus,j1*ones(size(aus)),'.','color',[0.5 0.5 0.5]);
end
ylabel('Cell ID')
% Distinguish assemblies
subplot(212);hold on;
for j1=1:nC
    aus = newSp(j1,:); aus(isnan(aus))=[];
    plot(aus,j1*ones(size(aus)),'.','color',[0.5 0.5 0.5]);
end
for j1=1:length(cellIDs)
    aus = spAs(j1,:); aus(isnan(aus))=[];
    plot(aus,cellIDs(j1)*ones(size(aus)),'.','MarkerSize',10);
end
ylabel('Cell ID')
xlabel('Time (s)')
%% Am I doing refractory period right?
min(diff(newSp,[],2),[],2)  % Should be <= tref
