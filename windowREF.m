function [refRIR,refH] = windowREF(Data,Tini,Tfin)
%[refRIR,refH] = windowREF(Data,Tini,Tfin) Applies Hanning window to the
%reference RIR between time given by Tini and Tfin. Obtains the windowed
%impulse response and the frequency response of the Reference Line.
%   Input:
%       - Data      : raw data. Structure
%       - Tini      : initial time. Scalar
%       - Tfin      : end time. Scalar
%   Output:
%       - refRIR    : windowed RIR. N x M
%       - refH      : windowed Frequency Response. Nf x M
%
% Author: Antonio Figueroa Durán
% Date: January 2023

%% ERROR HANDLING
% plotFlag default value
if nargin < 3, error('windowREF Error: Not enough input parameters.'), end
%% MAIN CODE
tHann = 1e-3;
NHann = Data.Fs*tHann/2;
h = hann(2*NHann);
hHalf1 = h(1:end/2);
hHalf2 = h(end/2+1:end);

T = [Tini Tfin];
N = floor(Data.Fs*T);
Nsamples = N(2)-N(1);

% Windowing - 1/2 Hanning (hann) window on each side
w = vertcat(zeros(N(1),1),hHalf1,...
    ones(Nsamples-2*NHann,1),hHalf2,...
    zeros(Data.Nsamples-Nsamples-N(1),1));

% Windowing Reference Line
refRIR = repmat(w,1,size(Data.Ref.h,2)).*Data.Ref.h;
[refH,~] = fftUniBi(refRIR);

end

