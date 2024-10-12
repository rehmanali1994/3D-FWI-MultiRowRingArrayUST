function [x_idx, y_idx, z_idx, ind] = ...
    sampled_cylinder(Nx, Ny, Nz, Nr, theta, Nrows, Nel)
%SAMPLED CYLINDER Creates a sample circle mask
%   Nx, Ny, Nz  -- Size of x, y, and z grids
%               -- (x, z) are in-plane coordinates
%               -- y is the elevational coordinate
%   Nr          -- Radius of circle in samples
%   theta       -- Sampled angles
%   Nrows       -- Number of rows of rings in cylinder
%   Nel         -- Samples between rows of rings in cylinder
    % In-plane coordinates
    [x_idx, z_idx, xz_ind] = sampled_circle(Nx, Nz, Nr, theta);
    % Out-of-plane coordinates
    y_idx_start = round((Ny-(Nel*(Nrows-1)+1))/2);
    y_idx = y_idx_start + Nel*(0:Nrows-1) + 1;
    % Combine coordinates
    [XZ_IND, Y_IDX] = meshgrid(xz_ind, y_idx);
    ind = XZ_IND + Nx*Nz*(Y_IDX-1);
end

function [x_idx, z_idx, ind] = sampled_circle(Nx, Nz, Nr, theta)
%SAMPLED CIRCLE Creates a sample circle mask
%   Nx, Nz -- Size of x and z grids
%   Nr -- Radius of circle in samples
%   theta -- Sampled angles
    x = (-(Nx-1)/2:(Nx-1)/2); z = (-(Nz-1)/2:(Nz-1)/2);
    x_idx = dsearchn(x(:), Nr*cos(theta(:)));
    z_idx = dsearchn(z(:), Nr*sin(theta(:)));
    ind = sub2ind([Nz, Nx], z_idx, x_idx);
end