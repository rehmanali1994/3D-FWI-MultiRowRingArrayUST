clear
clc

% Selector
fileID = 1; % 1 = left-breast MRI; 2 = right-breast MRI
viewID = 2; % 1 = transverse/cross-section; 2 = sagittal; 3 = coronal

% Select Simulation
switch fileID
    case 1
        filename = 'LeftBreastMRI.mat';
    case 2
        filename = 'RightBreastMRI.mat';
end

% Load 2D and 3D Reconstruction
filename2D = ['results2D/', filename];
load(filename2D);
filename3D = ['results3D/', filename];
load(filename3D);

% Flip Z-Axis for Breast MRI Visualization Purposes
VEL_ESTIM_2D_VOLUME = VEL_ESTIM_2D_VOLUME(:,:,end:-1:1); % 2D slice FWI
VEL_ESTIM = VEL_ESTIM(:,:,end:-1:1); % 3D FWI
VEL_TRUE = VEL_TRUE(:,:,end:-1:1); % Ground Truth

% Upsample 2D Slices to 3D
c_bkgnd = 1480;
[Xi, Yi, Zi] = meshgrid(xi, yi, zi);
VEL_ESTIM2D = interp3(xi, yi, zelem, ...
    single(VEL_ESTIM_2D_VOLUME), Xi, Yi, Zi, 'linear', c_bkgnd);

% Remove PML Before Display
PML = 15; % PML size to trim off
VEL_TRUE = VEL_TRUE(PML+1:end-PML,PML+1:end-PML,PML+1:end-PML);
VEL_ESTIM = VEL_ESTIM(PML+1:end-PML,PML+1:end-PML,PML+1:end-PML);
VEL_ESTIM2D = VEL_ESTIM2D(PML+1:end-PML,PML+1:end-PML,PML+1:end-PML);

% Visualize 2D and 3D Reconstructions Alongside Ground Truth
switch viewID
    case 1
        sliceViewer(permute(cat(2,VEL_TRUE,VEL_ESTIM2D,VEL_ESTIM),[1,2,3]));
    case 2
        sliceViewer(permute(cat(2,VEL_TRUE,VEL_ESTIM2D,VEL_ESTIM),[3,2,1]));
    case 3
        sliceViewer(permute(cat(1,VEL_TRUE,VEL_ESTIM2D,VEL_ESTIM),[3,1,2]));
end

% Adjust Display Settings
crange = [1440, 1580];
clim(crange);

% Calculate RMS error and Pearson's Correlation Coefficient
Ri = sqrt(Xi.^2 + Yi.^2); radius = 0.10;
Ri = Ri(PML+1:end-PML,PML+1:end-PML,PML+1:end-PML);
RMSE_3D = sqrt(mean((VEL_ESTIM(Ri<radius)-VEL_TRUE(Ri<radius)).^2));
RMSE_2D = sqrt(mean((VEL_ESTIM2D(Ri<radius)-VEL_TRUE(Ri<radius)).^2));
PCC_3D = corrcoef([VEL_TRUE(Ri<radius), VEL_ESTIM(Ri<radius)]);
PCC_3D = (PCC_3D(1,2)+PCC_3D(2,1))/2;
PCC_2D = corrcoef([VEL_TRUE(Ri<radius), VEL_ESTIM2D(Ri<radius)]);
PCC_2D = (PCC_2D(1,2)+PCC_2D(2,1))/2;