function windowedAudio = windowAudio(audioIn,Fs,varargin) 

% This function windows an audio file into shorter segments. It
% first truncates long silent periods, and then searches for short silent
% periods as natural break points within a pre-defined window length with
% tolerance. Finally, ramping is applied to the silent edges of each window. 

% Required arguments: audioIn, Fs

% Optional name-value pair arguments with defaults: 
% 'WindowLength',4        Length of the desired window in seconds
% 'SilenceThreshold',0.1  Threshold for silence detection, 0 < > 1
% 'RampLength',30      Length of ramping of silent edges in ms
% 'LongSilence',500       Minimum length of long silence to truncate in ms
% 'ShortSilence',100      Minimum length of short silence for split points in ms

defaultWin = 4;
defaultThresh = 0.1;
defaultRamp = 30;
defaultLong = 500;
defaultShort = 100;

p = inputParser;
addRequired(p,'audioIn',@isnumeric);
addRequired(p,'Fs',@isnumeric);
addParameter(p,'WindowLength',defaultWin,@checkWinLen);
addParameter(p,'SilenceThreshold',defaultThresh,@checkThresh);
addParameter(p,'RampLength',defaultRamp,@checkRamp);
addParameter(p,'LongSilence',defaultLong,@checkLong);
addParameter(p,'ShortSilence',defaultShort,@checkShort);

parse(p,audioIn,Fs,varargin{:});
audioIn  = p.Results.audioIn;
Fs = p.Results.Fs;
winLen  = p.Results.WindowLength;
silenceThreshold = p.Results.SilenceThreshold;
rampLength = p.Results.RampLength;
longSilence = p.Results.LongSilence/1000;
shortSilence = p.Results.ShortSilence/1000;

minWin = winLen-(0.5*winLen);
maxWin = winLen+(0.5*winLen);
msBuffer = .05 * Fs;

[n,m] = size(audioIn);
if n > 1 && m == 1
    %
elseif m > 1 & n == 1
    audioIn = audioIn';
elseif m == 2 && n > 2
    audioIn = mean(audioIn,2);
elseif n == 2 && m > 2
    audioIn = mean(audioIn,2)';
else
    error('Input audio should be a numeric vector containing 1 or 2 channels.')
end

%% Step 1. Cut long silences

maxLenSilence = Fs * longSilence;
env = imdilate(abs(audioIn), true(1501, 1));
silence = env < silenceThreshold;
silence = silence(:).';
[~,idx] = find(double(silence));

if idx(1) ~= 1
    silence(1:idx) = 1;
end
if idx(end) ~= 1
    silence(idx(end):end) = 1;
end
    
endSilence = [strfind(silence, [1 0]) numel(env)];
beginSilence = [1 strfind(silence, [0 1])];
silence = endSilence-beginSilence;
beginSilence(silence < maxLenSilence) = [];
endSilence(silence < maxLenSilence) = [];

beginSilence(2:end) = beginSilence(2:end)+msBuffer;
endSilence(1:end-1) = endSilence(1:end-1)-msBuffer;

toCut = [];
toCut(1:2:2*numel(beginSilence)) = beginSilence;
toCut(2:2:2*numel(beginSilence)+1) = endSilence;

dataOut = [audioIn zeros(numel(audioIn),1)];
        
for ii = 1:2:numel(toCut)
	dataOut(toCut(ii):toCut(ii+1),2) = 1;
end
dataOut(:,dataOut(end,:)==1) = [];
    
%% Step 2. Split windows at shorter silent points

maxLenSilence = Fs * shortSilence;
env = imdilate(abs(dataOut(:,1)), true(1501, 1));
silence = env < silenceThreshold;
silence = silence(:).';
[~,idx] = find(double(silence));
    
if idx(1) ~= 1
	silence(1:idx) = 1;
end
if idx(end) ~= 1
	silence(idx(end):end) = 1;
end

beginSilence = [1 strfind(silence, [0 1])];
endSilence = [strfind(silence, [1 0]) numel(env)];
silence = endSilence-beginSilence;
beginSilence(silence < maxLenSilence) = [];
endSilence(silence < maxLenSilence) = [];
beginSilence(2:end) = beginSilence(2:end)+msBuffer;
endSilence(1:end-1) = endSilence(1:end-1)-msBuffer;
    
toSplit = (beginSilence+endSilence)/2;
theWins = 1:(winLen*Fs):size(dataOut,1);
[~,~,~,theWins,~] = matchPoints(toSplit,theWins);
theWins = [1 round(theWins) size(dataOut,1)];

%% Step 3. Combine windows that are too short

idx = diff(theWins)<(minWin*Fs);
idx = find(logical([0 idx]));
for ii = 1:numel(idx)
    
    if idx(ii) == 2
        theWins(idx(ii)) = inf;
        
    elseif idx(ii)>2 && idx(ii) ~= numel(theWins)
        [~,appendTo] = min(diff([theWins(idx(ii)-1) theWins(idx(ii)) theWins(idx(ii)+1)]));
        if appendTo == 1
            theWins(idx(ii)) = inf;
        else
            theWins(idx(ii)+1) = inf;
        end
        
    elseif idx(ii) == numel(theWins)
        theWins(idx(ii)-1) = inf;
    end
end

theWins = theWins(~isinf(theWins));

%% Step 4. Split windows that are too long

idx = diff(theWins)>(maxWin*Fs);
idx = find(logical([0 idx]));
for ii = 1:numel(idx)
    newSplit = (theWins(idx(ii))-theWins(idx(ii)-1))/2;
    theWins(end+1) = theWins(idx(ii)-1) + newSplit;
end
theWins = round(sort(theWins));
theWins(theWins==size(dataOut,1)) = [];

%% Step 5. Compile

ramp = round(rampLength/1000*Fs);
rampVector = linspace(0,1,ramp)';

windowedAudio = cell(numel(theWins),1);
for ii = 1:numel(theWins)
    if ii ~= numel(theWins)
        winOut = dataOut(theWins(ii):theWins(ii+1),1);
    else
        winOut = dataOut(theWins(ii):end,1);
    end

    winOut(1:ramp) = winOut(1:ramp).*rampVector;
    winOut(end-ramp+1:end) = winOut(end-ramp+1:end).*rampVector(end:-1:1);
    
    windowedAudio{ii} = winOut;
    
end
end

%% Functions
function checkWinLen(x)
    if ~isnumeric(x)
        error('Window length should be a numeric input greater than 0 and less than the total duration (in seconds) of the audio data.');
    end
end
function checkThresh(x)
    if ~(x > 0 && x < 1)
        error('Silence threshold should be a numeric input greater than 0 and lesser than 1.');
    end
end
function checkRamp(x)
    if ~isnumeric(x)
        error('Ramping should be a numeric input in milliseconds.');
    end
end
function checkLong(x)
    if ~isnumeric(x)
        error('Long silence should be a numeric input in milliseconds.');
    end
end
function checkShort(x)
    if ~isnumeric(x)
        error('Short silence should be a numeric input in milliseconds.');
    end
end
function [eucDists,unmatched1,unmatched2,matched1,matched2] = matchPoints(vector1,vector2)

if numel(vector1)>0 && numel(vector2)>0
    
    if size(vector1,1) == 1
        vector1 = vector1';
    end
    
    if size(vector2,1) == 1
        vector2 = vector2';
    end
    
    if numel(vector1) > numel(vector2) % Longer vector becomes vectorB
        vectorA = vector2;
        vectorB = vector1;
        
        longerVector = 1;
        
    elseif numel(vector1) < numel(vector2)
        vectorB = vector2;
        vectorA = vector1;
        
        longerVector = 2;
        
    else
        vectorA = vector1;
        vectorB = vector2;
        
        longerVector = 0;
    end
    
    % Pair up observations
    counter = 1;
    eucDists = [];
    numObs = numel(vectorB);
    
    for i = 1:numObs
        if sum(~isinf(vectorA))~=0 && sum(~isinf(vectorB))~=0
            
            maxLen = max(numel(vectorB),numel(vectorA));
            padding = inf(maxLen,1);
            vectorA = [vectorA ; padding(1:abs(numel(vectorA)-maxLen))];
            vectorB = [vectorB ; padding(1:abs(numel(vectorB)-maxLen))];
            
            [closestA,idxA] = min(abs(vectorA(i)-vectorB));
            [closestB,idxB] = min(abs(vectorB(i)-vectorA));
            [~,ind] = max([idxA idxB]);
            
            if idxA == idxB % vectorA(i)/vectorB(i) are each other's closest match
                
                eucDists(counter) = closestA;
                
                hits(1,counter) = vectorA(i);
                hits(2,counter) = vectorB(i);
                
                vectorA(i) = inf;
                vectorB(i) = inf;
                
            else % there are potentially unpaired points
                switch ind
                    case 1
                        lower = vectorB(vectorB<vectorB(idxA));
                        if numel(lower) > 1
                            vectorB(vectorB<lower(end)) = inf;
                        end
                    case 2
                        lower = vectorA(vectorA<vectorA(idxB));
                        if numel(lower) > 1
                            vectorA(vectorA<lower(end)) = inf;
                            
                        end
                end
                
                compA = vectorA(~isinf(vectorA));
                compB = vectorB(~isinf(vectorB));
                maxElement = min(numel(compA),numel(compB));
                
                switch ind
                    case 1
                        [~,diffAB] = min(abs(compA(1)-compB));
                    case 2
                        [~,diffAB] = min(abs(compB(1)-compA));
                end
                
                diffAB = diffAB - 1;
                
                if maxElement == 1
                    compDiff = 1;
                elseif maxElement < 8 + diffAB
                    compA = compA(1:maxElement);
                    compB = compB(1:maxElement);
                    compDiff = maxElement-diffAB;
                else
                    compA = compA(1:8+diffAB);
                    compB = compB(1:8+diffAB);
                    compDiff = 8;
                end
                
                switch ind
                    
                    case 1 % vectorB(i) may be unpaired
                        
                        if sum(abs(compA(1:compDiff)-compB(1:compDiff))) ...
                                > sum(abs(compA(1:compDiff) ...
                                - compB(1+diffAB:compDiff+diffAB)))
                            
                            eucDists(counter) = closestA;
                            
                            hits(1,counter) = vectorA(i);
                            hits(2,counter) = vectorB(idxA);
                            
                            vectorA(i) = inf;
                            vectorB(idxA) = inf;
                            vectorB(i:idxA-1) = [];
                            
                        else
                            
                            if vectorB(i) ~= inf
                                
                                eucDists(counter) = abs(vectorA(i)-vectorB(i));
                                
                                hits(1,counter) = vectorA(i);
                                hits(2,counter) = vectorB(i);
                                
                                vectorA(i) = inf;
                                vectorB(i) = inf;
                                
                            else
                                % vectorB(i) is inf now!
                                counter = counter - 1;
                                vectorA = [vectorA(1:i) ; vectorA(i:end)];
                                vectorA(i) = inf;
                            end
                            
                        end
                        
                    case 2 % vectorA(i) may be unpaired
                        
                        if sum(abs(compB(1:compDiff)-compA(1:compDiff))) ...
                                > sum(abs(compB(1:compDiff) ...
                                - compA(1+diffAB:compDiff+diffAB)))
                            
                            eucDists(counter) = closestB;
                            
                            hits(1,counter) = vectorA(idxB);
                            hits(2,counter) = vectorB(i);
                            
                            vectorB(i) = inf;
                            vectorA(idxB) = inf;
                            vectorA(i:idxB-1) = [];
                            
                        else
                            
                            if vectorA(i) ~= inf
                                eucDists(counter) = abs(vectorB(i)-vectorA(i));
                                
                                hits(1,counter) = vectorA(i);
                                hits(2,counter) = vectorB(i);
                                
                                vectorA(i) = inf;
                                vectorB(i) = inf;
                                
                            else
                                counter = counter - 1;
                                vectorB = [vectorB(1:i) ; vectorB(i:end)];
                                vectorB(i) = inf;
                                
                            end
                            
                        end
                        
                end
            end
            
        end
        counter = counter + 1;
    end
    
    if longerVector == 1
        
        matched2 = hits(1,:);
        
        if sum(ismember(vector2,hits(1,:)))~=numel(vector2)
            unmatched2 = vector2(~ismember(vector2,hits(1,:)));
            unmatched2 = unmatched2(~isinf(unmatched2));
            if isempty(unmatched2)
                unmatched2 = 0;
            end
        else
            unmatched2 = 0;
        end
        
        matched1 = hits(2,:);
        
        if sum(ismember(vector1,hits(2,:)))~=numel(vector1)
            unmatched1 = vector1(~ismember(vector1,hits(2,:)));
            unmatched1 = unmatched1(~isinf(unmatched1));            
            if isempty(unmatched1)
                unmatched1 = 0;
            end
        else
            unmatched1 = 0;
        end
        
    else
        
        matched2 = hits(2,:);
        
        if sum(ismember(vector2,hits(2,:)))~=numel(vector2)
            unmatched2 = vector2(~ismember(vector2,hits(2,:)));            
            unmatched2 = unmatched2(~isinf(unmatched2));            
            if isempty(unmatched2)
                unmatched2 = 0;
            end
        else
            unmatched2 = 0;
        end

        matched1 = hits(1,:);
        
        if sum(ismember(vector1,hits(1,:)))~=numel(vector1)
            unmatched1 = vector1(~ismember(vector1,hits(1,:)));
            unmatched1 = unmatched1(~isinf(unmatched1));
            if isempty(unmatched1)
                unmatched1 = 0;
            end
        else
            unmatched1 = 0;
        end
    end
    
else
    eucDists = NaN;
    unmatched1 = NaN;
    unmatched2 = NaN;
    matched1 = NaN; 
    matched2 = NaN;
end

end