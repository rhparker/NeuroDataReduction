function [ETs] = genNHPP_wBreaks(rvec,tvec,breaktimes)
% Generate spike train when you need to break it up into 
%   high/low rate 
%
% 
ETs = [];
dt = tvec(2)-tvec(1);

% ASSUME breaktimes includes beginning, end of the desired time window!
nWin  = numel(breaktimes)-1;

for j1=1:nWin
    btime = breaktimes(j1); etime = breaktimes(j1+1);
    bind  = find(tvec>btime,1);eind = find(tvec>etime,1);
    
    ETtemp  = genNHPP(@(x)eval_pw_const(x,rvec(bind:eind),tvec(bind:eind)-btime),...
            etime-btime,1,dt);
     
    ETs = [ETs ETtemp+btime];    
end

