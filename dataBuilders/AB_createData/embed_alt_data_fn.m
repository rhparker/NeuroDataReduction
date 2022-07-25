function newdata = embed_alt_data_fn(spM,n,T,unit,scheme,varargin)
% This function will take ISI data, and embed original and compress
% patterns into it. The activity will be precisely sequential
% Input - n = [n1 n2] ; n1: number of original activity, n2: number of
% compress activity
% T=[T1 T2] - interspike interval of original and compress pattern
% unit - which neurons will be embeded, and their orders

% S = spM(unit,:);
nneu = size(spM,1);
assemblysize = numel(unit);
T1 = T(1); T2 = T(2);

n1 = n(1); n2 = n(2);

%generate uniformly random location
loc = sort(min(spM(:)) + rand(n1+n2,1)*max(spM(:)));

% divide loc for alternating scheme
if strcmp(scheme,'alt')
  loc1 = loc([1:n1/2 (n1+n2)/2+1:(n1+n2)/2+n1/2]);
  loc2 = setdiff(loc,loc1);
elseif strcmp(scheme,'twoblock')
  loc1 = loc(1:n1);
  loc2 = loc(n1+1:n1+n2);
elseif strcmp(scheme,'random')

end

X1 = [];
X2 = [];

c = 0; %jittering

%sequentially forward part
for i = 1:n1
  t = loc1(i)+ c*rand(1);
  y = zeros(assemblysize,1);
  for j = 1:assemblysize
    y(j) = t + (j-1)*T1; % a 0.2s apart sequential
  end
  X1 = [X1, y];
end

%forward and compressed
for i = 1:n2
  t = loc2(i)+ c*rand(1);
  y = zeros(assemblysize,1);
  for j = 1:assemblysize
    y(j) = t + (j-1)*T2; % a 0.04s apart sequential
  end
  X2 = [X2, y];
end
%% How to embed the X into spM data
%cell 46-50
S = spM(unit,:);
m = min(sum(~isnan(S),2));
rmloc = randperm(m,sum(n));
S(:,rmloc)=nan; %remove n spike

% Trim_Idx = find(sum(~isnan(S)),1,'last');

S(:,end-sum(n)+1:end) = [X1 X2];
S = sort(S,2);

% OtherUnit = setdiff([1:nneu],unit);

spM(unit,:) = S;


newdata.spM = spM;
newdata.acttime1 = X1;
newdata.acttime2 = X2;
newdata.assembly_ISI = T;
newdata.n = n;
newdata.unit = unit;




