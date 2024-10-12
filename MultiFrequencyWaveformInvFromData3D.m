clear
clc

% Load Functions
addpath(genpath('solvers4FWI/'));

% Load Simulation Information Created by GenKWaveSimInfo.m
save_filename = 'RightBreastMRI.mat';
siminfo_filename = ['sim_info/', save_filename];
load(siminfo_filename);

% Load Saved Dataset
dataset_filename = ['datasets/', save_filename];
load(dataset_filename, 'full_dataset', 'fDATA');

% Geometry
c_geom = 1500; % Sound Speed [mm/us] for Geometric TOF
zelem = (-(numRows-1)/2:(numRows-1)/2)*rowSpacing;
theta = -pi:2*pi/numElemPerRow:pi-2*pi/numElemPerRow;
[THETA, ZELEM] = meshgrid(theta, zelem);
transducerPositionsX = circle_radius*cos(THETA); 
transducerPositionsY = circle_radius*sin(THETA); 
transducerPositionsZ = ZELEM;

% Attenuation Iterations Always Happen After SoS Iterations
niterPerFreq = 5*ones(size(fDATA)); % Number of Sound Speed Iterations Per Frequency

% Geometry of Ground-Truth Simulation
xsim = (-(Nxi-1)/2:(Nxi-1)/2)*dxi;
ysim = (-(Nzi-1)/2:(Nzi-1)/2)*dzi;
zsim = (-(Nyi-1)/2:(Nyi-1)/2)*dyi;
[Xsim, Ysim, Zsim] = meshgrid(xsim, ysim, zsim);

% Window Time Traces - Extract the Frequencies for Waveform Inversion
REC_DATA = permute(full_dataset,[4,3,2,1]);
    % Tx Elem, Rx Row, Rx Elem, Frequencies

%% Waveform Inversion

% Create Sound Speed Map and Transducer Ring
dxi = 0.8e-3; xmax = 125.6e-3; 
dzi = 0.8e-3; zmax = 50.8e-3;
xi = -xmax:dxi:xmax; yi = xi; zi = -zmax:dzi:zmax; 
Nxi = numel(xi); Nyi = numel(yi); Nzi = numel(zi);
[Xi, Yi, Zi] = meshgrid(xi, yi, zi);
x_circ = transducerPositionsX;
y_circ = transducerPositionsY;
z_circ = transducerPositionsZ;
x_idx = dsearchn(xi(:), x_circ(:));
y_idx = dsearchn(yi(:), y_circ(:));
z_idx = dsearchn(zi(:), z_circ(:));
ind = sub2ind([Nyi, Nxi, Nzi], y_idx, x_idx, z_idx);
msk = zeros(Nyi, Nxi, Nzi); msk(ind) = 1;
x_idx = reshape(x_idx, [numRows, numElemPerRow]);
y_idx = reshape(y_idx, [numRows, numElemPerRow]);
z_idx = reshape(z_idx, [numRows, numElemPerRow]);
ind = reshape(ind, [numRows, numElemPerRow]);

% Ground Truth Sound Speed on Same Grid
VEL_TRUE = interp3(Xsim, Ysim, Zsim, ...
    medium.sound_speed, ...
    Xi, Yi, Zi, 'spline', c_geom);

% Parameters
alphaDB = 0.0; % Attenuation [dB/(MHz*mm)]
alphaNp = (log(10)/20)*alphaDB*((1e3)/(1e6)); % Attenuation [Np/(Hz*m)]
% Solve Options for Helmholtz Equation
sign_conv = -1; % Sign Convention
a0 = 20; % PML Constant
L_PML = 12.5e-3; % Thickness of PML  

% Conversion of Units for Attenuation Map
Np2dB = 20/log(10);
slow2atten = (1e6)/(1e2); % Hz to MHz; m to cm

% Which Subset of Transmits to Use
dwnsmp = 1; 
tx_include = 1:dwnsmp:numElemPerRow;
REC_DATA = REC_DATA(tx_include,:,:,:); 
REC_DATA(isnan(REC_DATA)) = 0; % Eliminate Blank Channel
NUM_BATCHES = 8; % HOW TO SPLIT THE BORN SERIES CALCULATION

% Extract Subset of Times within Acceptance Angle
numElemLeftRightExcl = 31; 
elemLeftRightExcl = -numElemLeftRightExcl:numElemLeftRightExcl;
elemInclude = true(numElemPerRow, numElemPerRow);
for tx_element = 1:numElemPerRow 
    elemLeftRightExclCurrent = elemLeftRightExcl + tx_element;
    elemLeftRightExclCurrent(elemLeftRightExclCurrent<1) = numElemPerRow + ...
         elemLeftRightExclCurrent(elemLeftRightExclCurrent<1);
    elemLeftRightExclCurrent(elemLeftRightExclCurrent>numElemPerRow) = ...
        elemLeftRightExclCurrent(elemLeftRightExclCurrent>numElemPerRow) - numElemPerRow;
    elemInclude(tx_element,elemLeftRightExclCurrent) = false;
end

% Remove Outliers from Observed Signals Prior to Waveform Inversion
perc_outliers = 0.99; % Confidence Interval Cutoff
for rx_row = 1:numRows
for f_idx = 1:numel(fDATA)
    REC_DATA_SINGLE_FREQ = squeeze(REC_DATA(:,rx_row,:,f_idx)); 
    signalMagnitudes = elemInclude(tx_include,:).*abs(REC_DATA_SINGLE_FREQ);
    num_outliers = ceil((1-perc_outliers)*numel(signalMagnitudes));
    [~,idx_outliers] = maxk(signalMagnitudes(:),num_outliers);
    REC_DATA_SINGLE_FREQ(idx_outliers) = 0;
    REC_DATA(:,rx_row,:,f_idx) = REC_DATA_SINGLE_FREQ;
end
end

% Initial Constant Sound Speed Map [m/s]
c_init = c_geom; % Initial Homogeneous Sound Speed [m/s] Guess
VEL_INIT = c_init*ones(Nyi,Nxi,Nzi); 
% Initial Constant Attenuation [Np/(Hz m)]
ATTEN_INIT = 0*alphaNp*ones(Nyi,Nxi,Nzi);

% Generate Sources
SRC = zeros(Nyi, Nxi, Nzi, numel(tx_include), 'single');
for elmt_idx = 1:numel(tx_include)
    % Single Element Source
    x_idx_src = x_idx(:,tx_include(elmt_idx)); 
    y_idx_src = y_idx(:,tx_include(elmt_idx)); 
    z_idx_src = z_idx(1,tx_include(elmt_idx)):z_idx(end,tx_include(elmt_idx)); 
    SRC(y_idx_src, x_idx_src, z_idx_src, elmt_idx) = 1; 
end

% (Nonlinear) Conjugate Gradient
search_dir = zeros(Nyi,Nxi,Nzi); % Conjugate Gradient Direction
gradient_img_prev = zeros(Nyi,Nxi,Nzi); % Previous Gradient Image
VEL_ESTIM = VEL_INIT; ATTEN_ESTIM = ATTEN_INIT;
SLOW_ESTIM = 1./VEL_ESTIM + ...
    1i*sign(sign_conv)*ATTEN_ESTIM/(2*pi); % Initial Slowness Image [s/m]
crange = [1440, 1580]; % For reconstruction display [m/s]
attenrange = 10*[-1,1]; % For reconstruction display [dB/(cm MHz)]
c0 = mean(VEL_ESTIM(:)); cutoff = 0.75; ord = Inf; % Parameters for Ringing Removal Filter
% Iterations of Conjugate Gradient
for f_idx = 1:numel(fDATA)
    % Iterations at Each Frequency
    for iter_f_idx = 1:niterPerFreq(f_idx)
        tic;
        % Step 1: Accumulate Backprojection Over Each Element
        iter = iter_f_idx + sum(niterPerFreq(1:f_idx-1));
        % Reset CG at Each Frequency (SoS and Attenuation)
        if (iter_f_idx == 1)
            search_dir = zeros(Nyi,Nxi,Nzi); % Conjugate Gradient Direction
            gradient_img_prev = zeros(Nyi,Nxi,Nzi); % Previous Gradient Image
        end
        % Forward Solve Helmholtz Equation
        [WVFIELD, VIRT_SRC] = solveHelmholtzBornSeries3D_Batches(xi, yi, zi, ...
            VEL_ESTIM, ATTEN_ESTIM, SRC, ...
            fDATA(f_idx), sign_conv, a0, L_PML, false, NUM_BATCHES);
        % Build Adjoint Sources
        scaling = zeros(numel(tx_include), 1);
        ADJ_SRC = zeros(Nyi, Nxi, Nzi, numel(tx_include));
        for elmt_idx = 1:numel(tx_include)
            WVFIELD_elmt = WVFIELD(:,:,:,elmt_idx);
            REC_SIM = WVFIELD_elmt(ind(:,elemInclude(tx_include(elmt_idx),:))); 
            REC = REC_DATA(elmt_idx, :, ...
                elemInclude(tx_include(elmt_idx),:), f_idx); 
            REC = squeeze(REC);
            scaling(elmt_idx) = (REC_SIM(:)'*REC(:)) / ...
                (REC_SIM(:)'*REC_SIM(:)); % Source Scaling
            ADJ_SRC_elmt = zeros(Nyi, Nxi, Nzi);
            ADJ_SRC_elmt(ind(:,elemInclude(tx_include(elmt_idx),:))) = ...
                scaling(elmt_idx)*REC_SIM - REC;
            ADJ_SRC(:,:,:,elmt_idx) = ADJ_SRC_elmt;
        end
        % Backproject Error
        ADJ_WVFIELD = solveHelmholtzBornSeries3D_Batches(xi, yi, zi, ...
            VEL_ESTIM, ATTEN_ESTIM, ADJ_SRC, ...
            fDATA(f_idx), sign_conv, a0, L_PML, true, NUM_BATCHES);
        SCALING = repmat(reshape(scaling, ...
            [1,1,1,numel(scaling)]), [Nyi, Nxi, Nzi, 1]);
        BACKPROJ = -real(conj(SCALING.*VIRT_SRC).*ADJ_WVFIELD);
        % Accumulate Gradient Over Each Transmitter
        gradient_img = sum(BACKPROJ,4);
        % Remove Ringing from Gradient Image
        gradient_img = ringingRemovalFilt3D(xi, yi, zi, ...
            gradient_img, c0, fDATA(f_idx), cutoff, ord);
        % Step 2: Compute New Conjugate Gradient Search Direction from Gradient
        % Conjugate Gradient Direction Scaling Factor for Updates
        if ((iter_f_idx == 1) || (iter_f_idx == 1+niterPerFreq(f_idx)))
            beta = 0; 
        else 
            betaPR = (gradient_img(:)'*...
                (gradient_img(:)-gradient_img_prev(:))) / ...
                (gradient_img_prev(:)'*gradient_img_prev(:));
            betaFR = (gradient_img(:)'*gradient_img(:)) / ...
                (gradient_img_prev(:)'*gradient_img_prev(:));
            beta = min(max(betaPR,0),betaFR);
        end
        search_dir = beta*search_dir-gradient_img;
        gradient_img_prev = gradient_img;
        % Step 3: Compute Forward Projection of Current Search Direction
        PERTURBED_WVFIELD = solveHelmholtzBornSeries3D_Batches(xi, yi, zi, ...
            VEL_ESTIM, ATTEN_ESTIM, VIRT_SRC.*search_dir, ...
            fDATA(f_idx), sign_conv, a0, L_PML, false, NUM_BATCHES);
        dREC_SIM = zeros(numRows, numElemPerRow, numel(tx_include));
        for elmt_idx = 1:numel(tx_include)
            % Forward Projection of Search Direction Image
            PERTURBED_WVFIELD_elmt = PERTURBED_WVFIELD(:,:,:,elmt_idx);
            dREC_SIM(:,elemInclude(tx_include(elmt_idx),:),elmt_idx) = -scaling(elmt_idx) * ...
                PERTURBED_WVFIELD_elmt(ind(:,elemInclude(tx_include(elmt_idx),:)));
        end
        % Step 4: Perform a Linear Approximation of Exact Line Search
        perc_step_size = 1; % (<1/2) Introduced to Improve Compliance with Strong Wolfe Conditions 
        alpha = -(gradient_img(:)'*search_dir(:))/(dREC_SIM(:)'*dREC_SIM(:));
        SLOW_ESTIM = SLOW_ESTIM + perc_step_size * alpha * search_dir;
        VEL_ESTIM = 1./real(SLOW_ESTIM); % Wave Velocity Estimate [m/s]
        ATTEN_ESTIM = 2*pi*imag(SLOW_ESTIM)*sign(sign_conv);
        % Visualize Numerical Solution - Plot Final Reconstructions
        rx_row = 16;
        subplot(1,3,1); imagesc(xi,yi,VEL_ESTIM(:,:,z_idx(rx_row,1)),crange);
        title(['Estimated Wave Velocity ', num2str(iter)]); axis image;
        xlabel('Lateral [m]'); ylabel('Axial [m]'); colorbar; colormap gray;
        subplot(1,3,2); imagesc(xsim,ysim,VEL_TRUE(:,:,z_idx(rx_row,1)),crange);
        title('True Wave Velocity'); axis image;
        xlabel('Lateral [m]'); ylabel('Axial [m]'); colorbar; colormap gray;
        disp(['Iteration ', num2str(iter)]); subplot(1,3,3); toc;
    end
end

% Plot Final Reconstructions
rx_row = 16;
subplot(1,2,1); imagesc(xi,yi,VEL_ESTIM(:,:,z_idx(rx_row,1)),crange);
title(['Estimated Wave Velocity ', num2str(iter)]); axis image;
xlabel('Lateral [m]'); ylabel('Axial [m]'); colorbar; colormap gray;
subplot(1,2,2); imagesc(xsim,ysim,VEL_TRUE(:,:,z_idx(rx_row,1)),crange);
title('True Wave Velocity'); axis image;
xlabel('Lateral [m]'); ylabel('Axial [m]'); colorbar; colormap gray;

% Save Final Image Reconstruction
save(['results/results3D/', save_filename], ...
    'xi', 'yi', 'zi', 'VEL_ESTIM', 'VEL_TRUE');