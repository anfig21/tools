function [b,theta,phi] = das(P,r,f,c,res,plotFlag)
%[b,theta,phi] = das(P,r,f,c,res,plotFlag) Obtain the Delay-and-Sum
%beamformer output for the microphone array located at r at the single
%frequency f.
%   Input:
%       - P         : frequency response at the specific frequency bins. M x 1
%       - r         : microphone positions. M x 3
%       - f         : frequency. Scalar
%       - c         : speed of sound in the medium. Scalar
%       - res       : angular resolution. Scalar
%       - plotFlag  : 'true' to plot beamformer power and max value (DOA)
%                     'false' (Default value)
%   Output:
%       - b         : beamformer output. Ntheta x Nphi
%       - theta     : elevation angles from zenith in rad. Ntheta x 1
%       - phi       : azimuth angles from X-axis in rad. Nphi x 1
%
% Author: Antonio Figueroa Durán
% Date: June 2023

%% ERROR HANDLING
if nargin < 5, error('das Error: Not enough input parameters.'), end
if nargin < 6, plotFlag = false; end

%% MAIN CODE
resRad = deg2rad(res);

% Steering coordinates
Ntheta = ceil(pi/resRad);   theta = linspace(0,pi,Ntheta);
Nphi = ceil(2*pi/resRad);   phi = linspace(0,(1-1/Nphi)*2*pi,Nphi);

% Steering vector
[Utheta,Uphi] = meshgrid(theta,phi);
UkSph = [ones(size(Utheta(:),1),1) Utheta(:) Uphi(:)];
[UkCart(:,1), UkCart(:,2), UkCart(:,3)] = sph2cart2(UkSph(:,1),UkSph(:,2),UkSph(:,3));

% Propagation matrix
K = 2*pi*f/c * UkCart;

% Delay vector
V = exp(1j*r*K');

% Beamformer output
b = reshape(V'*P,size(Uphi)).';

%% PLOT
% Beamformer output and max value
if plotFlag
    maximum = max(abs(b).^2,[],'all');
    [idxTheta,idxPhi] = find(abs(b).^2 == maximum);

    figure, hold on
    s = pcolor(rad2deg(phi),rad2deg(theta),abs(b).^2);
    plot(rad2deg(phi(idxPhi)),rad2deg(theta(idxTheta)),'Marker','x','MarkerSize',15)
    set(s,'edgecolor','none')
    xlabel('Aximuth $\phi$ / deg'), ylabel('Elevation $\Theta$ / deg')
    DR = colormapDR(abs(b).^2);
    colorbarpwn(0,DR)
    c = colorbar;
    applyColorbarProperties(c,'Beamformer power $|b(\theta,\phi)|^2$')
    applyAxisProperties(gca)
end

end