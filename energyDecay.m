function [eDecay,F0] = energyDecay(h,Fc,Fs,smoothing)
%[eDecay,F0] = energyDecay(h,Fs,smoothing) Obtains the energy decay curves
%of the set of impulse responses contained in h in third-octave bands.
%   Input:
%       - h         : impulse responses. nT x M
%       - Fc        : cut-off frequencies. 2 x 1
%       - Fs        : sampling rate. Scalar
%       - smoothing : smoothes the RIRs using moving average filter:
%                       'true'
%                       'false' (Default value)
%   Output:
%       - eDecay    : energy decay curves. nT x M x nB
%       - F0        : centre frequencies. numF0 x 1
%
% Author: Antonio Figueroa Durán
% Date: March 2023

%% ERROR HANDLING
if nargin < 3, error('energyDecay Error: Not enough input parameters.'), end
if nargin < 4 || isempty(smoothing), smoothing = false; end

%% MAIN CODE
[nT,M] = size(h);

% Third-octave band filter
Filter = octaveFilterBank('1/3 octave',Fs,'FrequencyRange',Fc,'FilterOrder',6);
F0 = getCenterFrequencies(Filter);
numF0 = length(F0);

% Removes the early part of the RIR
tEarly = 25e-3;
h(1:ceil(tEarly*Fs),:) = 0;

% Filtering
hf = Filter(h);
hf = permute(hf,[1 3 2]);

%% DECAY ESTIMATION
% Hilbert transform: gets the envelope
env = abs(hilbert(hf));

% Optional Moving Average filtering
if smoothing
    Nsamples = 101;
    B = ones(Nsamples,1)/Nsamples;
    movEnv = filtfilt(B,1,env);
    env = movEnv;
end

% Schroeder's Backward Integration
tD = 28300;     % Noise floor estimation
eDecay = zeros([nT M numF0]);
eDecay(1:tD,:,:) = flip(cumsum(flip(env(1:tD,:,:)).^2)./sum(env(1:tD,:,:).^2));

%% PLOT
% t = (1:nT)/Fs;
% ii = [1 6 9];   % Freq bands
% 
% figure, hold on
% for jj = 1:length(ii)
%     plot(t*1e3,10*log10(eDecay(:,1,ii(jj))))
%     plot(t*1e3,20*log10(env(:,1,ii(jj))/max(env(:,1,ii(jj)))))
% end
% xlabel('Time / ms'), ylabel('Energy decay')
% legend('$F_c=158$\,Hz','','$F_c=500$\,Hz','','$F_c=1$\,kHz','')
% applyAxisProperties(gca)
% applyLegendProperties(gcf)
% grid on

end

