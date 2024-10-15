function [wvfield, virtSrcs] = solveHelmholtzBornSeries3D(x, y, z, ...
    vel, atten, src, f, signConvention, a0, L_PML, adjoint)
    
% Discretization and Data Dimensions
dx = single(mean(diff(x))); 
dy = single(mean(diff(y))); 
dz = single(mean(diff(z))); 
Nx = numel(x); Ny = numel(y); Nz = numel(z); % Grid Points in X and Y
Nsrcs = size(src,4);

% Whether to Solve Adjoint Helmholtz Equation
if adjoint
    adjointSign = -1; % Change sign if adjoint
else
    adjointSign = 1; % Keep sign same if not adjoint
end

% Calculate Complex Slowness
SI = atten/(2*pi); SR = 1./vel;
S = SR + 1i*SI*sign(signConvention*adjointSign);
    
% Calculate Complex Wavenumber
k = 2*pi*f*S; % Wavenumber [1/m]

% Convert to gpuArray if GPU available
if canUseGPU
    k = gpuArray(k);
    src = single(gpuArray(src));
end

% Fourier grid - 2*pi*k actually
kx = 2*pi*(mod((0:Nx-1)/(dx*Nx)+1/(2*dx),1/dx)-1/(2*dx)); 
kx = reshape(kx,[1,Nx,1]);
ky = 2*pi*(mod((0:Ny-1)/(dy*Ny)+1/(2*dy),1/dy)-1/(2*dy)); 
ky = reshape(ky,[Ny,1,1]);
kz = 2*pi*(mod((0:Nz-1)/(dz*Nz)+1/(2*dz),1/dz)-1/(2*dz)); 
kz = reshape(kz,[1,1,Nz]);

% Definition of PML based on Born series paper
N = 9; % polynomial order
c = a0/L_PML; % attenuation [Np/m] inside PML
k0 = sqrt(mean(k.^2,'all')); % mean wavenumber [1/m]
f_boundary_curve = @(r) ...
    (c^2)*(N-c*r+2i*k0*r*signConvention*adjointSign).*((c*r).^(N-1)) ./ ...
    (factorial(N)*polyval(1./factorial(N:-1:0),c*r));
% Apply PML to obtain modified wavenumber map
x_pml = subplus(abs(x)+L_PML-(Nx-1)*dx/2);
x_pml = reshape(x_pml,[1,Nx,1]);
y_pml = subplus(abs(y)+L_PML-(Ny-1)*dy/2);
y_pml = reshape(y_pml,[Ny,1,1]);
z_pml = subplus(abs(z)+L_PML-(Nz-1)*dz/2);
z_pml = reshape(z_pml,[1,1,Nz]);
PML = f_boundary_curve(sqrt(x_pml.^2+y_pml.^2+z_pml.^2));
k = sqrt(k.^2 + PML);


% Construct potential map V: 
%   [ V(r) = k(r).^2 - k_0.^2 - i*epsilon ]
% where epsilon must satisfy:
%   [ epsilon >= max(|k(r).^2 - k_0.^2|)  ]
k_0 = (min(real(k(:)))+max(real(k(:))))/2; % medium wavenumber 
V = k.^2-k_0.^2; % Scattering Potential
epsilon = max(abs(V(:)))*signConvention*adjointSign;
V = V-1i*epsilon; % Full potential map V
gamm = 1i/epsilon*V; % Preconditioner gamma

% How many iterations needed for propagate over domain
%   [ pseudo-propagation length per iteration = 2*k_0/epsilon ]
pseudoproplen = 2*k_0/abs(epsilon); % pseudo-propagation length
max_distance = sqrt((Ny*dy)^2+(Nx*dx)^2+(Nz*dz)^2); % max distance over grid
max_iterations = ceil(max_distance/pseudoproplen); 
fprintf('Maximum iterations = %d\n', max_iterations);

% 3D FFT/IFFT - Plan to replace with 3D cuFFT in future
fft3 = @(sig) fft(fft(fft(sig,[],3),[],2),[],1); % fftn
ifft3 = @(sig) ifft(ifft(ifft(sig,[],3),[],2),[],1); % ifftn

% Calculate Greens function (in Fourier domain/k-space)
%   [ g_0(p) = 1/(|p|.^2 - k_0^2 - i*epsilon) ]
p2 = kx.^2 + ky.^2 + kz.^2; % k-space [|p|.^2]
g0_k = 1./(p2-(k_0^2 + 1i*epsilon)); % This one can be precalculated
G = @(csrc) ifft3(g0_k .* fft3(csrc)); % Green's function operator

% Allocate memory 
if canUseGPU
    wvfield = zeros(Ny,Nx,Nz,Nsrcs,'single','gpuArray'); % Final Solution
else
    wvfield = zeros(Ny,Nx,Nz,Nsrcs,'single'); % Final Solution
end

% Start simulation loop
updateBornSeries = @(wvf) wvf-gamm.*(wvf-G(V.*wvf-src)); % Born-Series Update
for it = 1:max_iterations
    % Convergent Born Series Solution 
    wvfield = updateBornSeries(wvfield); % Update wavefield
%    % Show Propagation
%    elmt_idx = 1;
%    slice_z_idx = 64; % which element's wavefield to show
%    Eshow = wvfield(:,:,slice_z_idx,elmt_idx);
%    Eshow = real(squeeze(Eshow));
%    imagesc(x,y,Eshow,[-1,1]*(1e-15));
%    title(['wavefiled solution Iter = ', num2str(it)]);
%    axis square; colorbar;
%    xlabel('x (\lambda)','FontSize',16);
%    ylabel('y (\lambda)','FontSize',16);
%    set(gca,'FontSize',14); drawnow;
end

% Compute [dH/ds u] where H is Helmholtz matrix and u is the wavefield
if nargout > 1 % Linear Indexing Into Sparse Array
    % Create Sparse Matrix
    wvfield_matrix = (4*pi*f)*k;
    % Compute Virtual Sources
    virtSrcs = wvfield_matrix.*wvfield;
    % Move from GPU to CPU
    if canUseGPU
        virtSrcs = gather(virtSrcs);
    end
end

% Move from GPU to CPU
if canUseGPU
    wvfield = gather(wvfield);
end

end
