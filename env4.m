function envelopeOut = env4(signalIn,FsIn,FsOut)

% Tilsen, S., & Arvaniti, A. (2013). Speech rhythm analysis with
% decomposition of the amplitude envelope: characterizing rhythmic patterns
% within and across languages. The Journal of the Acoustical Society of
% America, 134(1), 628-639.

% Written by Sam Tilsen.

par.bandpass = [400 4000];  %: vocalic energy bandpass range (there is no a priori justification for this range)
par.lowpass = 10;           %: lowpass filtering cutoff
par.Fs = FsIn;                %: sampling rate

%extract the envelope:
env = envm_band_energy(signalIn,par);

%normalize:
env = env-mean(env);
env = env/max(abs(env));

% downsample
fsRatio = round(FsIn/FsOut);
envelopeOut = resample(env,1,fsRatio);

end
function be = envm_band_energy(wav,par)

%INPUTS
% wav: the acoustic signal
% par.bandpass: vocalic energy bandpass range 
% par.lowpass: lowpass filtering cutoff
% par.Fs: sampling rate
% par.ds: downsampling factor

%OUTPUTS
%band energy (e.g. vocalic energy  amplitude envelope), unnormalized

if ~isfield(par,'bandpass_order'), par.bandpass_order = 4; end
if ~isfield(par,'lowpass_order'), par.lowpass_order = 4; end

%PARAMETERS
Fs              = par.Fs;
bandpass        = par.bandpass;                           %bandpass frequencies
lowpass         = par.lowpass;                            %lowpass cutoff frequency
nyquist         = Fs/2;                                     %nyquist frequency
bandpass_norm   = bandpass/nyquist;                         %normalized bandpass vector
bandpass_order  = par.bandpass_order;
lowpass_norm    = lowpass/nyquist;                          %normalized lowpass cutoff
lowpass_order   = par.lowpass_order;

%PROCESSING
wav_dc          = wav - mean(wav);                          %remove DC
[bbp,abp]       = butter(bandpass_order,bandpass_norm);     %1st order, bandpass filter pareters
sig_bp          = filtfilt(bbp,abp,wav_dc);                   %bandpass filter
%sig_power      = (abs(sig_bp).^2);                         %power
sig_mag         = (abs(sig_bp));                            %magnitude
[blp,alp]       = butter(lowpass_order,lowpass_norm);       %4th order, lowpass (<20Hz) filter pareters
sig_lp          = filtfilt(blp,alp,sig_mag);                  %lowpass filter
sig_dc          = (sig_lp - mean(sig_lp));                  %remove DC

be = sig_dc;

end

