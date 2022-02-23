function envelopeOut = env1(signalIn,FsIn,FsOut,lowPass,winSz)

% Jarne, C. (2017). A heuristic approach to obtain signal envelope with a
% simple software implementation. arXiv preprint arXiv:1703.06812.

% This method uses an arbitrary window size to determine the moving
% maxima. The default to this optional argument is 250 samples.

% Set window sample size
if exist('winSz','var')
    if isempty(winSz)
        winSz = 250;
    end
elseif ~exist('winSz','var')
    winSz = 250;
end

env = abs(signalIn);
numWin = length(env)/winSz;
if mod(numWin,1) ~= 0
    numWin = numWin - 1;
end

% Moving max
for i = 1:winSz:numel(env)
    idx = i:i+winSz-1;
    if i == numWin
        idx = i:numel(env);
    end
    idx(idx>numel(env)) = [];
    env(idx) = max(env(idx));
end

% Low pass filter
[z,p,k] = butter(10,lowPass/(round(FsIn/2)),'low');
[sos,g] = zp2sos(z,p,k);
env = filtfilt(sos,g,env);

% Resample
[p,q] = rat(FsOut/FsIn); 
env = resample(env,p,q);

envelopeOut = env;
end