function SpikeCount = pair_Generator(T,p,m,L,binw,indtype)
% Generate a pairt of spike counts with T/binw bins
% type = 'indep', 'sync'
p_low = p(1);
p_high = p(2);
SpikeCount = nan(2,T/binw);
T_bin = nan(1,T/binw);
max_int = floor(T/L)-1;

if m >max_int
  err('m is too large \n');
end

% select random m intervals
int_I = datasample([1:1:max_int],m,'Replace',false);

% flag the elementary bin with high prob
T_Index = zeros(T,1);
for i = int_I
  T_Index(i*L:(i+1)*L-1) = 1;
end

% General Spike Train
Sp1 = binornd(1,p_low,T,1).*(~T_Index) + binornd(1,p_high,T,1).*T_Index;

if strcmp(indtype,'sync')
  Sp2 = binornd(1,p_low,T,1).*(~T_Index) + binornd(1,p_high,T,1).*T_Index;
elseif strcmp(indtype,'indep')
  % select random m intervals
  int_I = datasample([1:1:max_int],m,'Replace',false);
  % flag the elementary bin with high prob
  T_Index = zeros(T,1);
  for i = int_I
    T_Index(i*L:(i+1)*L-1) = 1;
  end
  Sp2 = binornd(1,p_low,T,1).*(~T_Index) + binornd(1,p_high,T,1).*T_Index;
end

% Convert to Spike Count
for j = 1: T/binw
  T_bin(1,j) = sum(T_Index((j-1)*binw+1:j*binw));
  SpikeCount(1,j) = sum(Sp1((j-1)*binw+1:j*binw));
  SpikeCount(2,j) = sum(Sp2((j-1)*binw+1:j*binw));
end
% %%
% subplot(2,2,1)
% plot(T_bin);
% subplot(2,2,2)
% plot(SpikeCount(1,:));
% subplot(2,2,3)
% plot(SpikeCount(2,:));
clearvars Sp1 Sp2 T_Index;