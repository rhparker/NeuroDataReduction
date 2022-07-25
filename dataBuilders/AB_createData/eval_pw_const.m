function [f] = eval_pw_const(t,fvec,tvec)
% eval_pw_const: evaluate a piecewise constant function
%
%  tvec: vector of time values
%  fvec: function values
%  Rule: If tvec(j)<=t<tvec(j+1), then f = 

if numel(t)>1
    f = zeros(size(t));
    for j1=1:numel(t)
        whichj = find(t(j1)<tvec,1);

        if (isempty(whichj))
            % t is larger than any value
            f(j1) = fvec(end);
        else
            f(j1) = fvec(whichj);
        end
    end
else
    whichj = find(t<tvec,1);

    if (isempty(whichj))
        % t is larger than any value
        f = fvec(end);
    else
        f = fvec(whichj);
    end
end
    
end

