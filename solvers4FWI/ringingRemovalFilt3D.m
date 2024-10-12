function img_out = ringingRemovalFilt3D(xin, yin, zin, img_in, c0, f, cutoff, ord)
%RINGINGREMOVALFILT Filter to Remove Ringing from Waveform Inversion Result
%
% img_out = ringingRemovalFilt(xin, yin, zin, img_in, c0, f, cutoff)
% INPUT:
%   xin = x (N-element array) grid [in m] for the input image
%   yin = y (M-element array) grid [in m] for the input image
%   zin = z (K-element array) grid [in m] for the input image
%   img_in = input image (M x N x K array)
%   c0 = reference sound speed [m/s]
%   f = frequency [Hz]
%   cutoff = ranging from 0 (only keep DC) to 1 (up to ringing frequency)
%   ord = order of radial Butterworth filter cutoff (Inf for sharp cutoff)
% OUTPUT:
%   img_out = output image (M x N x K array)

% Number and spacing of input points (assumes a uniform spacing)
Nxin = numel(xin); dxin = mean(diff(xin));
Nyin = numel(yin); dyin = mean(diff(yin));
Nzin = numel(zin); dzin = mean(diff(zin));

% K-Space
kxin = fftshift(((0:Nxin-1)/Nxin)/dxin);
kxin(kxin>=1/(2*dxin)) = kxin(kxin>=1/(2*dxin)) - 1/dxin;
kyin = fftshift(((0:Nyin-1)/Nyin)/dyin);
kyin(kyin>=1/(2*dyin)) = kyin(kyin>=1/(2*dyin)) - 1/dyin;
kzin = fftshift(((0:Nzin-1)/Nzin)/dzin);
kzin(kzin>=1/(2*dzin)) = kzin(kzin>=1/(2*dzin)) - 1/dzin;
[Kxin, Kyin, Kzin] = meshgrid(kxin, kyin, kzin);

% Cutoff Wavenumber
k = 2*f/c0; kcutoff = k*cutoff;

% Filter Based on Cutoff
radialFilt = 1./(1+((Kxin.^2 + Kyin.^2 + Kzin.^2)/(kcutoff^2)).^ord);

% Apply Filter
img_out = real(ifftn(ifftshift(radialFilt.*(fftshift(fftn(img_in))))));

end