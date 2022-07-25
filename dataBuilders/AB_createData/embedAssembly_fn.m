function [new_spM] = embedAssembly_fn(spAs,spM,cellID)
%
% spAs      = assembly spikes
% spM       = background spikes
% cellID    = cell IDs for assmebly

% Refractory period
tref = 0.015;

% Check validity of arguments
nC      = size(spM,1);
nCas    = size(spAs,1);

if (max(cellID)>nC || min(cellID<=0))
    error('cellID not valid: either negative, or not enough background cells');
end

% Check for repeats?
if (~all(diff(sort(cellID))>0))
    error('cellID not valid: there is a repeated index');
end

% Maximal possible size
new_spM = nan(nC,size(spM,2)+size(spAs,2));

for j1=1:nCas
    % Background spikes
    temp    = spM(cellID(j1),:);
    % Assembly spikes
    tempAs  = spAs(j1,:);
    
    tempAs(isnan(tempAs))=[];  %Remove nan
    for j2=1:length(tempAs)
        % For each spike, remove any background spikes within tref
        % Use "nan" for this
        z = tempAs(j2);
        temp(temp>z-tref & temp<z+tref) =nan;
    end
    % Add new spikes
    temp = [temp tempAs];
    
    % And now remove nan
    temp = sort(temp);
    
    new_spM(cellID(j1),1:length(temp))=temp;  

end
% For unchanged cells, copy over remaining 
uncin = setdiff(1:nC,cellID);
for j1=uncin
    new_spM(j1,1:size(spM,2)) = spM(j1,:);
end
  
% Delete unneccessary columns
new_spM(:, all(isnan(new_spM)))=[];

end

