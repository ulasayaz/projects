function [xdata] = normalize_frame(xdata)
% converts to double and normalizes a matrix
% to double and normalise to 1
xdata=double(xdata);
xdata=xdata-min(xdata(:));
xdata=xdata/max(xdata(:));


end

