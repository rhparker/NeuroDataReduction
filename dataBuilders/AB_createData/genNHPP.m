
function [ EventTimes ] = genNHPP(rate_fh,T,n,dt)
%GENNHPP Generate a NHPP
%   Generates a NHPP, N(t): # of events by time t
% INPUTS
%   rate_fh : function handle for rate (vectorized)
%   T : end of time horizon
%   n : number of sample paths to generate (default n = 1)
%   dt : expected time-scale of rate variations
%
% OUTPUTS
%   EventTimes : if n = 1, EventTimes is a vector (length is number of events)
%                if n > 1, EventTimes is a cell array with n rows (# columns is number of events)
%   NumArrived : Number Arrived by time t (time is row, columns are sample paths)

% Input Error Checking ****************************************************
narginchk(2, 4)
if nargin < 3 || isempty(n), n = 1; end % Default is 1 sample path
if nargin < 4 || isempty(dt), dt = 1; end % Default is 1 day
if ~isa(rate_fh,'function_handle')
    error('rate_fh must be a valid function handle')
end
if T <= 0, error('T must be a positive real number'), end
% End (Input Error Checking) **********************************************

% Generate a NHPP 
MaxLambda = max(rate_fh([0:dt:T])); 
 % ^^^^^This can be a very expensive call, depending on how
 % rate_fh is implemented!!!!!
 %
 
 % To cut overhead, do division in place
 %
%ph=@(t) rate_fh(t)/MaxLambda;

%Preallocate
EventTimes = nan(1,round(T*MaxLambda));
if n ==1 
    % Single Sample Path
    t = 0; Nevents = 0; done = false;
    
    while ~done
        t = t + (-1/MaxLambda)*log(rand());
        if t > T
            done = true;
        else
            if (rand()<= (rate_fh(t)/MaxLambda))  %ph(t);
                Nevents = Nevents+1;
                EventTimes(Nevents) = t;
            end
        end
    end
    % Cut off "nan" part 
    EventTimes(isnan(EventTimes))=[];
else 
    % Multiple Sample Paths
    NumPaths = n; EventTimes = cell(NumPaths,1);
    for path = 1:NumPaths
        tEvents = nan(1,round(T*MaxLambda));  %Preallocate
        t = 0; Nevents = 0; done = false; 
        while ~done
            t = t + (-1/MaxLambda)*log(rand());
            if t > T
                done = true;
            else
                if rand() <= ph(t);
                    Nevents = Nevents+1;
                    tEvents(Nevents)=t;
                    %EventTimes{path,Nevents} = t;
                end
            end
        end
        % Cut off "nan" part 
        tEvents(isnan(tEvents))=[];
        EventTimes{path} = tEvents;
    end
end



