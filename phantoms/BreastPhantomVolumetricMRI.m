function [C, c_bkgnd] = BreastPhantomVolumetricMRI(xq, yq, zq, leftBreast)
%BREASTPHANTOMVOLUMETRICMRI Outputs Sound Speed Phantom
% [C, c_bkgnd] = BreastPhantomVolumetricMRI(xq, yq, zq, leftBreast)
% INPUTS:
%   xq, yq, zq  -- grid over which sound speed is defined
%               -- (x,z) are in-plane coordinates
%               -- y is elevation dimension
%   leftBreast  -- boolean (use left [true] or right [false] breast)
% OUTPUTS:
%   C           -- sound speed map on grid [m/s]
%   c_bkgnd     -- background sound speed [m/s]
%   rho         -- density map on grid [kg/m^3]
%   d_bkgnd     -- background density [kg/m^3]

% Load 3D Skull Volume
load('BreastPhantomFromMRI.mat', 'breast_mri');
breast_mri = double(breast_mri);

% Coordinates
tmp = yq; yq = zq; zq = tmp; % swap y and z coordinates
xq = 1000*xq; yq = 1000*yq; zq = 1000*zq; % [m] to [mm]
nxq = numel(xq); nyq = numel(yq); nzq = numel(zq); 

% Pixel Size
xLen = 340; yLen = 340; zLen = 158.4; % Lengths of X,Y,Z Dimensions [mm]
dx = xLen/size(breast_mri,1); % set dx from image--in plane [mm]
dy = yLen/size(breast_mri,2); % set dy from image--in plane [mm]
dz = zLen/size(breast_mri,3); % MRI slice thickness [mm]

% Physical Grid [mm] for 3D Breast Volume
x = dx*(-(size(breast_mri,1)-1)/2:(size(breast_mri,1)-1)/2);
y = dy*(-(size(breast_mri,2)-1)/2:(size(breast_mri,2)-1)/2); 
z = dz*(-(size(breast_mri,3)-1)/2:(size(breast_mri,3)-1)/2);
[X, Y, Z] = meshgrid(x,y,z); 

% Indices for Breast Center on Grid
if leftBreast
    xc_point = 330; yc_point = 125; zc_point = 75;
else
    xc_point = 120; yc_point = 130; zc_point = 70;
end
% Skull Center in Physical Coordinates on Grid [mm]
xc = x(xc_point); yc = y(yc_point); zc = z(zc_point);

% Coordinates of Box to Interpolate Over Before Rotation
[Xq_unrot, Yq_unrot, Zq_unrot] = meshgrid(xq,yq,zq);
Vq_unrot = [Xq_unrot(:)'; Yq_unrot(:)'; Zq_unrot(:)'];

% Rotate the Box We are Interpolating Over 
if leftBreast
    theta = 100; % Angle [deg] from x 
    phi = 90; % Angle [deg] from z
else
    theta = 80; % Angle [deg] from x 
    phi = 80; % Angle [deg] from z
end
rot_matrix = [cosd(theta)*cosd(phi), -sind(theta), cosd(theta)*sind(phi); ... 
              sind(theta)*cosd(phi),  cosd(theta), sind(theta)*sind(phi); ...
                         -sind(phi),            0,             cosd(phi)];
Vq_rot = rot_matrix * Vq_unrot;

% Interpolate 3D Breast Volume Over Our Rotated Box
shrinkFactor = 1;
Xq = shrinkFactor*reshape(Vq_rot(1,:),[nxq, nyq, nzq]);
Yq = shrinkFactor*reshape(Vq_rot(2,:),[nxq, nyq, nzq]);
Zq = shrinkFactor*reshape(Vq_rot(3,:),[nxq, nyq, nzq]); 
extrap_val = 0; % Extrapolate Into Air in MRI Image Units
breast_seg = interp3(X, Y, Z, breast_mri, ...
    Xq+xc, Yq+yc, Zq+zc, 'cubic', extrap_val);

% Segmentation: Convert Air Outside the Breast into Water 
thresh = 40; % threshold estimate based on automatic threshold
tissue_pix = zeros(size(breast_seg),'logical');
for i = 1:size(breast_seg,3) % generate binary maps of the thresholds
    tissue_pix(:,:,i) = imfill(imbinarize(breast_seg(:,:,i), thresh),'holes'); 
end 

% Convert to Sound Speed Map
if leftBreast
    c_min = 1400; % Minimum Sound Speed [m/s]
    c_max = 1750; % Maximum Sound Speed [m/s]
else
    c_min = 1350; % Minimum Sound Speed [m/s]
    c_max = 1750; % Maximum Sound Speed [m/s]
end
c_bkgnd = 1500; % Sound Speed [m/s] in Water Outside the Breast
breast_radius_max = 80; % Outside this Breast Radius [mm] is water
C = (c_max-c_min)*(breast_seg-min(breast_seg(tissue_pix))) / ...
    (max(breast_seg(tissue_pix))-min(breast_seg(tissue_pix))) + c_min;
C(~tissue_pix) = c_bkgnd;
C(Xq_unrot.^2+Yq_unrot.^2>breast_radius_max^2) = c_bkgnd; 

end