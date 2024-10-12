clear
clc

% Load Functions
addpath(genpath('../k-Wave'));
addpath(genpath(pwd));

% Sound Speed Map 
gridSpacing = 0.4e-3; % grid spacing (dx, dy, dz) in each dimension [m]
dxi = gridSpacing; dzi = gridSpacing; dyi = gridSpacing;
PMLSize = 10; % Size of PML in pixels
PMLAlpha = 1000; % Absorption of PML [Np/m]
Nxi = 556; Nzi = 556; Nyi = 196; % Decided based on PMLSize
    % Keep each dimension+2*PMLSize some product of 2s, 3s, 5s, and 7s
xi = (-(Nxi-1)/2:(Nxi-1)/2)*dxi; % in-plane dimension [m]
yi = (-(Nyi-1)/2:(Nyi-1)/2)*dyi; % elevation [m]
zi = (-(Nzi-1)/2:(Nzi-1)/2)*dzi; % in-plane dimension [m]
[Xi, Zi, Yi] = meshgrid(xi, zi, yi);
leftBreast = false; % which breast to use, left [true] or right [false]?
[C, c_bkgnd] = BreastPhantomVolumetricMRI(xi, yi, zi, leftBreast);

% Attenuation Map
atten_bkgnd = 0.0025; % Background Attenuation [dB/(MHz^y cm)]
sos2atten = 2.5e-3; % Conversion from Sound Speed Difference to Attenuation
atten_varying = sos2atten*abs(C-c_bkgnd); % Varying Attenuation [dB/(MHz^y cm)]
atten = (atten_bkgnd + atten_varying); % Total Attenuation [dB/(MHz^y cm)]
y_atten = 1.01; % Power Law of Attenuation [y]
    % y cannot exactly equal 1 without the following line:
    % medium.alpha_mode = 'no_dispersion'

% Create Multi-Row Transducer Ring
circle_radius = 110e-3; 
numElemPerRow = 256; % Number of Ring Elements Per Row
numRows = 32; % Number of Rows in Transducer Ring
rowSpacingPixels = 6; % Pixels Between Rows in Ring
rowSpacing = rowSpacingPixels*dyi; % Spacing Between Rows in Ring
circle_rad_pixels = floor(circle_radius/gridSpacing);
theta = -pi:2*pi/numElemPerRow:pi-2*pi/numElemPerRow;
x_circ = circle_radius*cos(theta); 
z_circ = circle_radius*sin(theta); 
[x_idx, y_idx, z_idx, ind] = sampled_cylinder(Nxi, Nyi, Nzi, ...
    circle_rad_pixels, theta, numRows, rowSpacingPixels);
%msk = zeros(Nzi, Nxi, Nyi); msk(ind) = 1;

% Define the Properties of the Propagation Medium
rho_bkgnd = 1000; % Density [kg/m^3]
impedance = c_bkgnd*rho_bkgnd; % Impedance
medium.sound_speed = C; % Sound Speed [m/s]
medium.density = (rho_bkgnd/c_bkgnd)*C; % Proportional Density [kg/m^3]
medium.alpha_coeff = atten; % [dB/(MHz^y cm)]
medium.alpha_power = y_atten; % cannot exactly equal 1 without:
                              % medium.alpha_mode = 'no_dispersion'
medium.alpha_mode = 'no_dispersion'; % IGNORE VELOCITY DISPERSION!
kgrid = kWaveGrid(Nzi, dzi, Nxi, dxi, Nyi, dyi); % K-Space Grid Object

% Create Time Array
t_end = sqrt((Nzi*dzi)^2+(Nxi*dxi)^2+(Nyi*dyi)^2)/double(min(C(:))); 
%t_end = sqrt((2*circle_radius)^2+(rowSpacing*numRows)^2)/double(min(C(:))); 
cfl = 0.3; [kgrid.t_array, dt] = makeTime(kgrid, double(C), cfl, t_end);
fs = 1/dt; % Sampling Frequency [Hz]

% Define Properties of the Tone Burst Used to Drive the Transducer
% Transmit Pulse and Impulse Response
fracBW = 0.99;  % Fractional Bandwidth of the Transducer
fTx = 0.6e6; % Transmit Frequency [Hz]
tc = gauspuls('cutoff', fTx, fracBW, -6, -80); % Cutoff time at -100dB, fracBW @ -6dB
t = (-ceil(tc/dt):1:ceil(tc/dt))*dt; % (s) Time Vector centered about t=0
src_amplitude = 1e1;
N = 5; [b,a] = butter(N,2*fTx*dt*(1+fracBW));
gauss_pulse = filtfilt(b,a,gauspuls(t,fTx,fracBW));
dgauss_pulse = (gauss_pulse([2:end,1])-gauss_pulse([end,1:end-1]))/(2*pi*fTx*dt);
d2gauss_pulse = (dgauss_pulse([2:end,1])-dgauss_pulse([end,1:end-1]))/(2*pi*fTx*dt);
tx_signal = src_amplitude * d2gauss_pulse; % Calculate Transmit Pulse
emission_samples = numel(tx_signal); 
t_offset = -(emission_samples-1)*dt/2;
% Time Axis and Fouier Transform
time = kgrid.t_array+t_offset; % Time Axis [s]
fmax = 1.5e6; % Maximum Frequency [Hz]
df = 4e3; % Frequency Increments [Hz] 
          % -- No more than 1/(duration of simulation [s])
freq = 0:df:fmax; % Freqeuncies in Waveform [Hz]

% Define a Binary Sensor Mask
sensor.mask = zeros(Nzi, Nxi, Nyi);
sensor.mask(ind) = 1; % Make Sensor the Same As Source

% Create a Display Mask to Display the Transducer
display_mask = sensor.mask;

% Assign the Input Options
input_args = {'DisplayMask', display_mask, 'PMLInside', false, ...
    'PlotPML', false, 'PMLAlpha', PMLAlpha, 'PMLSize', PMLSize...
    'PlotSim', false, 'DataCast', 'gpuArray-single'};

% Save Simulation Info to File
if leftBreast
    filename = 'sim_info/LeftBreastMRI.mat';
else
    filename = 'sim_info/RightBreastMRI.mat';
end
save(filename, '-v7.3', 'numElemPerRow', 'numRows', ...
    'Nxi', 'Nyi', 'Nzi', 'dxi', 'dyi', 'dzi', ...
    'circle_radius', 'rowSpacing', ...
    'x_idx', 'y_idx', 'z_idx', 'ind', 'tx_signal', ...
    'time', 'freq', 'kgrid', 'medium', 'sensor', 'input_args');