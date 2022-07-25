function spM = generating_ISI_fn(ncell,T,baserate)

%% parameters
% ncell = 10;
D = 0.9*eye(ncell);
v = 0.2;
% T = 2000; %second - for 10 datasets section
sigma = 0.01;
tref = 15/100;

rbar = baserate*ones(ncell,1);% base firing rate per seconds

s0 = zeros(ncell,1);
t = 0;
ISI = [];
counter = 0;

%% generating

while t<T
  counter = counter + 1;
  %%
  s1 = D*s0 + normrnd(zeros(ncell,1),sigma*ones(ncell,1));
  %rate
  r = (1 + erf(v*(s1 - mean(s1))./sigma)).*rbar;
  
  %draw the time from exponential distribution, and add a refractory period
  spt = exprnd(1/r);
  
  % add to the data matrix
  ISI = [ISI spt'];
  
  %check
  t = min(sum(ISI,2));
  
  s0 = s1;
  
end
% 
spM = nan(ncell,max(sum(ISI~=0,2)));

%place ISI into spM
for j = 1:ncell
  temp = ISI(j,:)~=0;
  spM(j,1:sum(temp)) = ISI(j,temp);
end


%% add the refactory period
spM = spM + tref;


%% collect spike stamps
spM = cumsum(spM,2);
spM(spM>T) = nan;

%cut out the nan
spM = spM(:,1:max(sum(~isnan(spM),2)));
