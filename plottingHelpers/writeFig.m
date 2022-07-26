function [] = writeFig(fhandle, ahandle,Tpath,ext,ftyp,res)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
r = 400;
if ~strcmp(ftyp, 'pdf')
    set(findall(fhandle,'-property','FontName'),'FontName','Helvetica');
    set(findall(ahandle,'-property','FontName'),'FontName','Helvetica');
end

%%% Nature Neuro wants Helvetica Fonts.  How does this work?
    set(findall(fhandle,'-property','FontName'),'FontName','Helvetica');
    set(findall(ahandle,'-property','FontName'),'FontName','Helvetica');

    
%%% Matlab fig
if ~exist( sprintf('%s/mFigs', Tpath), 'dir')
    mkdir( sprintf('%s/mFigs', Tpath) );
end
savefig(fhandle,...
        sprintf('%s/mFigs/%s.fig', Tpath, ext),...
        'compact');
    
%%% exported graphics
exportgraphics(ahandle, ...
               sprintf('%s/%s.%s', Tpath, ext, ftyp), ...
               'Resolution',res, 'ContentType', 'vector' )
end

