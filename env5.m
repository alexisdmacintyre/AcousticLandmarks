function envelopeOut = env5(signalIn,FsIn,FsOut,lowPass)

% This method splits the broadband signal via a gammatone filter bank,
% applies a power law to the filter output, and re-combines to form a
% single envelope.

% Biesmans, W., Das, N., Francart, T., & Bertrand, A. (2016).
% Auditory-inspired speech envelope extraction methods for improved
% EEG-based auditory attention detection in a cocktail party scenario. IEEE
% Transactions on Neural Systems and Rehabilitation Engineering, 25(5),
% 402-412.

%% Parameters 

if exist('FsOut','var')
    if isempty(FsOut)
        FsOut = 1000;
    end
elseif ~exist('FsOut','var')
    FsOut = 1000;
end

if exist('lowPass','var')
    if isempty(lowPass)
        lowPass = 10;
    end
elseif ~exist('lowPass','var')
    lowPass = 10;
end

[pRat,q] = rat(FsOut/FsIn); 

[z,p,k] = butter(8,lowPass/(round(FsIn/2)),'low');
[sos,g] = zp2sos(z,p,k);

gammaFiltBank = gammatoneFilterBank('SampleRate',FsIn, ...
                                    'NumFilters',28, ...
                                    'FrequencyRange',[50,5000]);
                                
%% Generate envelope                                

audioOut = gammaFiltBank(signalIn);

audioOut = abs(audioOut);
audioOut = audioOut.^0.6;

env = mean(audioOut,2);

env = filtfilt(sos,g,env);

envelopeOut = rescale(resample(env,pRat,q));

end