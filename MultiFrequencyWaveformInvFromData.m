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
transducerPositionsX = circle_radius*cos(theta); 
transducerPositionsY = circle_radius*sin(theta); 
numElements = numel(transducerPositionsX);
clearvars theta

% Attenuation Iterations Always Happen After SoS Iterations
niterPerFreq = 5*ones(size(fDATA)); % Number of Sound Speed Iterations Per Frequency

% Find Corresponding Slice of Ground-Truth
Nel = rowSpacing/dyi;
y_idx_start = round((Nyi-(Nel*(numRows-1)+1))/2);
y_idx = y_idx_start + Nel*(0:numRows-1) + 1;
VEL_TRUE = medium.sound_speed(:,:,y_idx);
xsim = (-(Nxi-1)/2:(Nxi-1)/2)*dxi;
ysim = (-(Nzi-1)/2:(Nzi-1)/2)*dzi;

%% Waveform Inversion Settings

% Create Sound Speed Map and Transducer Ring
dxi = 0.8e-3; xmax = 125.6e-3;
xi = -xmax:dxi:xmax; yi = xi;
Nxi = numel(xi); Nyi = numel(yi);
[Xi, Yi] = meshgrid(xi, yi);
x_circ = transducerPositionsX;
y_circ = transducerPositionsY;
x_idx = dsearchn(xi(:), x_circ(:));
y_idx = dsearchn(yi(:), y_circ(:));
ind = sub2ind([Nyi, Nxi], y_idx, x_idx);
msk = zeros(Nyi, Nxi); msk(ind) = 1;

% Parameters
alphaDB = 0.0; % Attenuation [dB/(MHz*mm)]
alphaNp = (log(10)/20)*alphaDB*((1e3)/(1e6)); % Attenuation [Np/(Hz*m)]
% Solve Options for Helmholtz Equation
sign_conv = -1; % Sign Convention
a0 = 20; % PML Constant
L_PML = 12.5e-3; % Thickness of PML  

% Extract Subset of Times within Acceptance Angle
numElemLeftRightExcl = 31; %31; %63; % 127;
elemLeftRightExcl = -numElemLeftRightExcl:numElemLeftRightExcl;
elemInclude = true(numElements, numElements);
for tx_element = 1:numElements 
    elemLeftRightExclCurrent = elemLeftRightExcl + tx_element;
    elemLeftRightExclCurrent(elemLeftRightExclCurrent<1) = numElements + ...
         elemLeftRightExclCurrent(elemLeftRightExclCurrent<1);
    elemLeftRightExclCurrent(elemLeftRightExclCurrent>numElements) = ...
        elemLeftRightExclCurrent(elemLeftRightExclCurrent>numElements) - numElements;
    elemInclude(tx_element,elemLeftRightExclCurrent) = false;
end

% Initial Constant Sound Speed Map [m/s]
c_init = c_geom; % Initial Homogeneous Sound Speed [m/s] Guess
VEL_INIT = c_init*ones(Nyi,Nxi); 
% Initial Constant Attenuation [Np/(Hz m)]
ATTEN_INIT = 0*alphaNp*ones(Nyi,Nxi);

% Which Subset of Transmits to Use
dwnsmp = 1; 
tx_include = 1:dwnsmp:numElements;

%% Loop of 2D FWI over each row of receive signals

% 2D FWI Over Each Slice in Volume
VEL_ESTIM_2D_VOLUME = zeros(Nyi, Nxi, numRows);

% Loop over each row of receive signals
for rx_row = 1:numRows

fsa_rf_data = squeeze(full_dataset(:,:,rx_row,:));
REC_DATA = permute(fsa_rf_data,[3,2,1]);
REC_DATA = REC_DATA(tx_include,:,:); 

% Remove Outliers from Observed Signals Prior to Waveform Inversion
perc_outliers = 0.99; % Confidence Interval Cutoff
for f_idx = 1:numel(fDATA)
    REC_DATA_SINGLE_FREQ = REC_DATA(:,:,f_idx); 
    signalMagnitudes = elemInclude(tx_include,:).*abs(REC_DATA_SINGLE_FREQ);
    num_outliers = ceil((1-perc_outliers)*numel(signalMagnitudes));
    [~,idx_outliers] = maxk(signalMagnitudes(:),num_outliers);
    REC_DATA_SINGLE_FREQ(idx_outliers) = 0;
    REC_DATA(:,:,f_idx) = REC_DATA_SINGLE_FREQ;
end

% (Nonlinear) Conjugate Gradient
search_dir = zeros(Nyi,Nxi); % Conjugate Gradient Direction
gradient_img_prev = zeros(Nyi,Nxi); % Previous Gradient Image
VEL_ESTIM = VEL_INIT; ATTEN_ESTIM = ATTEN_INIT;
SLOW_ESTIM = 1./VEL_ESTIM + ...
    1i*sign(sign_conv)*ATTEN_ESTIM/(2*pi); % Initial Slowness Image [s/m]
crange = [1350, 1600]; % For reconstruction display [m/s]
attenrange = 10*[-1,1]; % For reconstruction display [dB/(cm MHz)]
c0 = mean(VEL_ESTIM(:)); cutoff = 0.75; ord = Inf; % Parameters for Ringing Removal Filter
for f_idx = 1:numel(fDATA)
    % Iterations at Each Frequency
    for iter_f_idx = 1:niterPerFreq(f_idx)
        % Step 1: Accumulate Backprojection Over Each Element
        iter = iter_f_idx + sum(niterPerFreq(1:f_idx-1));
        % Reset CG at Each Frequency (SoS and Attenuation)
        if (iter_f_idx == 1)
            search_dir = zeros(Nyi, Nxi); % Conjugate Gradient Direction
            gradient_img_prev = zeros(Nyi, Nxi); % Previous Gradient Image
        end
        % Generate Sources
        SRC = zeros(Nyi, Nxi, numel(tx_include));
        for elmt_idx = 1:numel(tx_include)
            % Single Element Source
            x_idx_src = x_idx(tx_include(elmt_idx)); 
            y_idx_src = y_idx(tx_include(elmt_idx)); 
            SRC(y_idx_src, x_idx_src, elmt_idx) = 1; 
        end
        % Forward Solve Helmholtz Equation
        [WVFIELD, VIRT_SRC] = solveHelmholtzBornSeries(xi, yi, VEL_ESTIM, ATTEN_ESTIM, ...
            SRC, fDATA(f_idx), sign_conv, a0, L_PML, false);
        % Build Adjoint Sources
        scaling = zeros(numel(tx_include), 1);
        ADJ_SRC = zeros(Nyi, Nxi, numel(tx_include));
        for elmt_idx = 1:numel(tx_include)
            WVFIELD_elmt = WVFIELD(:,:,elmt_idx);
            REC_SIM = WVFIELD_elmt(ind(elemInclude(tx_include(elmt_idx),:))); 
            REC = REC_DATA(elmt_idx, ...
                elemInclude(tx_include(elmt_idx),:), f_idx); REC = REC(:);
            scaling(elmt_idx) = (REC_SIM(:)'*REC(:)) / ...
                (REC_SIM(:)'*REC_SIM(:)); % Source Scaling
            ADJ_SRC_elmt = zeros(Nyi, Nxi);
            ADJ_SRC_elmt(ind(elemInclude(tx_include(elmt_idx),:))) = ...
                scaling(elmt_idx)*REC_SIM - REC;
            ADJ_SRC(:,:,elmt_idx) = ADJ_SRC_elmt;
        end
        % Backproject Error
        ADJ_WVFIELD = solveHelmholtzBornSeries(xi, yi, VEL_ESTIM, ATTEN_ESTIM, ...
            ADJ_SRC, fDATA(f_idx), sign_conv, a0, L_PML, true);
        SCALING = repmat(reshape(scaling, [1,1,numel(scaling)]), [Nyi, Nxi, 1]);
        BACKPROJ = -real(conj(SCALING.*VIRT_SRC).*ADJ_WVFIELD);
        % Remove Ringing from Gradient Image
        gradient_img = ringingRemovalFilt(xi, yi, sum(BACKPROJ,3), c0, fDATA(f_idx), cutoff, ord);
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
        PERTURBED_WVFIELD = solveHelmholtzBornSeries(xi, yi, VEL_ESTIM, ATTEN_ESTIM, ...
            VIRT_SRC.*search_dir, fDATA(f_idx), sign_conv, a0, L_PML, false);
        dREC_SIM = zeros(numel(tx_include), numElements);
        for elmt_idx = 1:numel(tx_include)
            % Forward Projection of Search Direction Image
            PERTURBED_WVFIELD_elmt = PERTURBED_WVFIELD(:,:,elmt_idx);
            dREC_SIM(elmt_idx,elemInclude(tx_include(elmt_idx),:)) = -permute(scaling(elmt_idx) * ...
                PERTURBED_WVFIELD_elmt(ind(elemInclude(tx_include(elmt_idx),:))),[2,1]);
        end
        % Step 4: Perform a Linear Approximation of Exact Line Search
        perc_step_size = 1; % (<1/2) Introduced to Improve Compliance with Strong Wolfe Conditions 
        alpha = -(gradient_img(:)'*search_dir(:))/(dREC_SIM(:)'*dREC_SIM(:));
        SLOW_ESTIM = SLOW_ESTIM + perc_step_size * alpha * search_dir;
        VEL_ESTIM = 1./real(SLOW_ESTIM); % Wave Velocity Estimate [m/s]
        ATTEN_ESTIM = 2*pi*imag(SLOW_ESTIM)*sign(sign_conv);
        % Plot Final Reconstructions at this Slice
        subplot(1,2,1); imagesc(xi,yi,VEL_ESTIM,crange);
        title(['Estimated Wave Velocity Slice ', num2str(rx_row)]); axis image;
        xlabel('Lateral [m]'); ylabel('Axial [m]'); colorbar; colormap gray;
        subplot(1,2,2); imagesc(xsim,ysim,VEL_TRUE(:,:,rx_row),crange);
        title(['True Wave Velocity Slice ', num2str(rx_row)]); axis image;
        xlabel('Lateral [m]'); ylabel('Axial [m]'); 
        colorbar; colormap gray; drawnow;
    end
end

% Add Slice Into Reconstructed Volume
VEL_ESTIM_2D_VOLUME(:,:,rx_row) = VEL_ESTIM;

end

% Save Final Image Reconstruction
save(['results/results2D/', save_filename], 'xi', 'yi', 'zelem', 'VEL_ESTIM_2D_VOLUME', 'VEL_TRUE');