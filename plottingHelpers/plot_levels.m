function [] = plot_levels(xdata, ydata, newfig, varargin)
% =========================================================================
% Plots bars as boxes with heights given by ydata and box widths centered
% at xdata
% newfig = true if you want to create a new figure
% newfig = false if you want to add to existing axis
% varargin = all plotting options that could be passed to plot()
% =========================================================================


ydata = ydata(:)';
ydata_new = [ydata; ydata];
ydata_new = ydata_new(:);

xdata = xdata(:)';

xdata_diff = diff( xdata );

xdata_diff = [xdata_diff(1) xdata_diff xdata_diff(end)];

xdata_new = [xdata - xdata_diff(1:end-1)/2; ...
             xdata + xdata_diff(2:end)/2];

xdata_new = xdata_new(:);

if newfig
    figure
end

plot(xdata_new, ydata_new, varargin{:})
end

