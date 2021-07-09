function [lowX,midX] = findSignalCrossings(x,minThresh)

if exist('minThresh','var')
    if isempty(minThresh)
        minThresh = 0.005;
    end
elseif ~exist('minThresh','var')
    minThresh = 0.005;
end

% Find positive going signal
y = [diff(x); 0];
y(y<0) = 0;
y(y>0) = 1;
y([1 numel(y)]) = 0;
l = strfind(y', [0 1])';
u = strfind(y', [1 0])';

lowX = l+((u-l)/4);
midX = l+((u-l)/2);

% Threshold = amplitude gained over how many samples

len = u-l;
am = x(u)-x(l);
rt = am./len;
lowX(rt<minThresh) = [];
midX(rt<minThresh) = [];

end