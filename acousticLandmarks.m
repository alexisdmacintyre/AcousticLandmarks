clearvars
close all
clc

%% Acoustic Landmark Identification

% This script accompanies: 

% Pushing the envelope: evaluating speech rhythm with different envelope
% extraction techniques, MacIntyre, A.D., Cai, C.Q., and Scott, S.K.
% (2021).

% Written by Alexis Deighton MacIntyre.

% Additional Contributions:

% Code for Envelope 3 is adapted from Oganian, Y., & Chang, E. F. (2019). A speech
% envelope landmark for syllable encoding in human superior temporal gyrus.
% Science advances, 5(11).

% Code for Envelope 4 is adapted from Tilsen, S., & Arvaniti, A.
% (2013). Speech rhythm analysis with decomposition of the amplitude
% envelope: characterizing rhythmic patterns within and across languages.
% The Journal of the Acoustical Society of America, 134(1), 628-639.

%% Parameters

% Acoustic Landmark
spEnvelope = 3; % 1. Moving max. 2. Hilbert transform. 3. Bark scale. 4. Vocalic energy.
sigEvent = 5; % 1. Lower crossings. 2. Mid-crossings. 3. Peaks. 4. Bases of peaks. 5. Peaks in the first derivative.

% Preprocessing and Windowing
cleanUp = 1; % Set to 1 to preprocess WAV (this will probably improve the windowing process)
winSp = 1; % Set to 1 to perform windowing
winLen = 3; % Length of window in seconds
rampLen = 50; % Length of ramping to apply to window edges in ms

% Signal Events
eventProm = 0.1; % Prominence parameter for peak finding
eventDist = 80; % Inter-peak distance parameter for peak finding (in sample Fs)
eventThresh = 0.005; % Threshold parameter for lower and mid-crossing finding

%% Load WAV

localPath = fileparts(matlab.desktop.editor.getActiveFilename);
[speechWav,Fs] = audioread(fullfile(localPath,'example.wav'));

%% Preprocess audio

if cleanUp == 1
    CV = @(alpha) -sqrt(2) * erfcinv(2*alpha);
    alpha = 0.9995;
    zs = CV([(1-alpha)/2   1-(1-alpha)/2]);
    speechWav = detrend(speechWav); % Detrend
    speechWav = speechWav-mean(speechWav); % Remove offset
    speechWav = speechWav./(max(abs(speechWav))); % Normalise
    sigLim = mean(speechWav) + zs*std(speechWav); % Limit
    speechWav(speechWav<sigLim(1)) = sigLim(1);
    speechWav(speechWav>sigLim(2)) = sigLim(2);
    speechWav = rescale(speechWav,-1,1); % Rescale
end

%% Window audio

if winSp == 1
    speechOut = windowAudio(speechWav,Fs,'WindowLength',winLen,'RampLength',rampLen);
else
    speechOut = {speechWav};
end

%% Extract envelope

FsOut = 1000; % Envelope sample rate
lowPass = 10; % 10 Hz

for i = 1:numel(speechOut)
    switch spEnvelope
        case 1
            featureOut = env1(speechOut{i,1},Fs,FsOut,lowPass);
        case 2
            featureOut = env2(speechOut{i,1},Fs,FsOut,lowPass);
        case 3
            featureOut = env3(speechOut{i,1},Fs,FsOut); % Oganian and Chang, 2019
        case 4
            featureOut = env4(speechOut{i,1},Fs,FsOut); % Tilsen and Arvaniti, 2013
    end
        
    % Smooth noise in signal floor and remove spurious peaks
    levels = statelevels(featureOut,round(FsOut/5),'mean');
    levels(1) = levels(1)-(0.25*levels(1));
    featureOut(featureOut<=levels(1)) = levels(1);
    
    featureOut = rescale(featureOut,-1,1);
    speechOut{i,2} = featureOut;
end

%% Identify Signal Events

for i = 1:size(speechOut,1)
    
    featureIn = speechOut{i,2};
    
    % For lower and mid-crossings, split into bi-level signal
    levs = [quantile(featureIn,.49) quantile(featureIn,.51)];
    if levs(1) == levs(2)
        jj = 0.1;
        while levs(1) == levs(2)
            levs = [quantile(featureIn,.49) quantile(featureIn,.51+jj)];
            jj = jj + 0.1;
        end
    end
    
    switch sigEvent
        case 1 % Lower Crossings
            [eventsOut,~] = findSignalCrossings(featureIn,eventThresh);
            
        case 2 % Mid-crossings
            [~,eventsOut] = findSignalCrossings(featureIn,eventThresh);

        case 3 % Peaks
            peaks = islocalmax(featureIn,'MinProminence',eventProm,...
                'MinSeparation',eventDist,'FlatSelection','first');
            eventsOut = find(peaks);
            
        case 4 % Bases of peaks
            peaks = islocalmax(featureIn,'MinProminence',eventProm,...
                'MinSeparation',eventDist,'FlatSelection','first');
            peaks = find(peaks);
            eventsOut = findSlopeBases(featureIn,peaks);
            
        case 5 % Peaks in first derivative
            
            chRt = diff(featureIn);
            chRt(chRt<0) = 0;
            chRt = rescale(chRt,-1,1);
            
            peaks = islocalmax(chRt,'MinProminence',eventProm,...
                'MinSeparation',eventDist,'FlatSelection','first');
            peaks = find(peaks);
            eventsOut = peaks;
    end
    
    speechOut{i,3} = eventsOut;
    
end

%% Plot Results

id = 1; % indices of window to plot

spSig = rescale(speechOut{id,1},-1,1);
feSig = speechOut{id,2};
evs = speechOut{id,3};

figure('Position',[300 300 1000 800],'Color','w');

clear p
ax1 = subplot(2,1,1);
p(1) = plot(spSig,'Color',[.88 .88 .88],'LineWidth',1);

ax2 = subplot(2,1,2);
p(2) = plot(feSig,'Color','k','LineWidth',1.25);
hold on
p(3) = plot(evs,feSig(round(evs)),'d',...
            'Color','k','LineWidth',1,'MarkerSize',10);

[f,e] = makeLabels(spEnvelope,sigEvent);
l = legend([p(1) p(2) p(3)],...
    'Speech Wave Form',f,e,...
    'EdgeColor','none','FontSize',14,...
    'Location','northoutside');

ts = 1:(numel(spSig)/Fs)/10:numel(spSig)/Fs;
ts = round(ts,2);
tix1 = 1:(numel(spSig)/numel(ts)):numel(spSig);
tix2 = 1:(numel(feSig)/numel(ts)):numel(feSig);


set(ax1,'XTick',tix1,'YLim',[-1.1 1.1],'XLim',[1 numel(spSig)])  
set(ax2,'XTick',tix2,'YLim',[-1.1 1.1],'XLim',[1 numel(feSig)])  
set(ax1,'YColor','none','YTick',[],'YGrid','off')
set(ax2,'YColor','none','YTick',[],'YGrid','off')
set(ax1,'XTickLabel',ts,'FontSize',12)
set(ax2,'XTickLabel',ts,'FontSize',12)

xlabel('Time (s)', 'FontSize',14)

%% Functions

function [f,e] = makeLabels(spEnvelope,signalEvent)

if spEnvelope == 1
    f = 'Envelope (moving max)';
elseif spEnvelope == 2
    f = 'Envelope (Hilbert transform)';
elseif spEnvelope == 3
    f = 'Envelope (Bark scale)';
elseif spEnvelope == 4
    f = 'Envelope (vocalic energy)';
end
    
if signalEvent == 1
    e = 'Lower crossings';
elseif signalEvent == 2
    e = 'Mid-crossings';
elseif signalEvent == 3
    e = 'Peaks';
elseif signalEvent == 4
    e = 'Bases of peaks';
elseif signalEvent == 5
    e = 'Peaks in the first derivative';
end

end