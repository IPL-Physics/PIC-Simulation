function [a, phi_grid, E_grid] = getAcc(pos, Nx, boxsize, n0, Gmtx, Lmtx)
    % Calculate the acceleration on each particle due to electric field
    N = numel(pos);
    dx = boxsize / Nx;
    j = floor(pos / dx) + 1;
    jp1 = j + 1;
    weight_j = (jp1 * dx - pos) / dx;
    weight_jp1 = (pos - j * dx) / dx;
    jp1 = mod(jp1 - 1, Nx) + 1;   % periodic BC

    n = accumarray(j, weight_j, [Nx, 1]) + accumarray(jp1, weight_jp1, [Nx, 1]);
    n = n * n0 * boxsize / N / dx;

    % Solve Poisson's Equation: laplacian(phi) = n-n0
    phi_grid = Lmtx \ (n - n0);

    % Apply Derivative to get the Electric field
    E_grid = -Gmtx * phi_grid;

    % Interpolate grid value onto particle locations
    E = weight_j .* E_grid(j) + weight_jp1 .* E_grid(jp1);

    a = -E;
end