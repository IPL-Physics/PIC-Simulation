function plasma_pic_simulation(Na,Nb,Nc,tEnd,va,vb,vc,vtha,vthb,vthc)
    % Simulation parameters
    %Na,Nb, and Nc are the counts of three separate distributions of
    %electrons. tEnd is how long the simulation will go for. va,vb, and vc
    %are beam velocities or mean velocity of each of the separate
    %distributions. and vtha,vthb, and vthc are the variance or thermal
    %velocities of those distributions
    N = Na + Nb + Nc;
    Nx = 400;     % Number of mesh cells
    time = 0;       % current time of the simulation
    dt = 1;       % timestep
    boxsize = 50;      % periodic domain [0, boxsize]
    n = 1;
    n0 = 1;       % electron number density
    df = 4;
    A = 0.0;   % perturbation
    % number of timesteps
    Nt = ceil(tEnd / dt);
    storeValues = true; % switch on for storing velocity data over interval timesteps
    truncation = 1; % skip timesteps to collect velocity data for limited sample (1 means collecting all data)
    Vel = nan(N, ceil(Nt/truncation));
    Pos = nan(N, ceil(Nt/truncation));
    Efield = zeros(Nx,Nt);
    start = 1; % what timestep to start storing velocities
    
    % Plotting parameters
    plotElectricField = true;
    plotPos = true;
    plotRealTime = true;
    
    % Generate Initial Conditions
    rng(42); % set the random number generator seed
    % construct 2 opposite-moving Gaussian beams
    posa = rand(Na, 1) * boxsize;
    posb = rand(Nb, 1) * boxsize;
    posc = rand(Nc, 1) * boxsize;
    pos = [posa; posb; posc];
    vela = vtha * randn(Na, 1);
    velb = vthb * randn(Nb, 1) + vb;
    velc = vthc * randn(Nc, 1) + vc;
    vel = [vela; velb; velc];
    
    % add perturbation
    velb = velb .* (1 + A * sin(2 * pi * posb / boxsize));
    velc = velc .* (1 + A * sin(2 * pi * posc / boxsize));
    
    % Construct matrix G to compute Gradient (1st derivative)
    dx = boxsize / Nx;
    e = ones(Nx, 1);
    Gmtx = spdiags([-e e], [-1 1], Nx, Nx);
    Gmtx(1, Nx) = -1;
    Gmtx(Nx, 1) = 1;
    Gmtx = Gmtx / (2 * dx);
    
    % Construct matrix L to compute Laplacian (2nd derivative)
    Lmtx = spdiags([e -2*e e], [-1 0 1], Nx, Nx);
    Lmtx(1, Nx) = 1;
    Lmtx(Nx, 1) = 1;
    Lmtx = Lmtx / dx^2;
    
    % calculate initial electric accelerations
    [acc, phi_grid] = getAcc(pos, Nx, boxsize, n0, Gmtx, Lmtx);
    
    % prep figure
    figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);

for i = 1:Nt
    % (1/2) kick
    vel = vel + acc * dt * 0.5;

    % drift (and apply periodic boundary conditions)
    pos = pos + vel * dt;
    pos = mod(pos, boxsize);

    % Suppress the warning for singular matrices
    warning('off', 'MATLAB:nearlySingularMatrix');

    % update accelerations
    [acc, phi_grid, E_grid] = getAcc(pos, Nx, boxsize, n0, Gmtx, Lmtx);

    % Turn the warning back on
    warning('on', 'MATLAB:nearlySingularMatrix');

    % (1/2) kick
    vel = vel + acc * dt * 0.5;

    % update time
    time = time + dt;

    if storeValues
        if i >= start && mod(i, truncation) == 0
            Vel(:, floor(i / truncation) - floor(start / truncation) + 1) = vel;
            Pos(:, floor(i / truncation) - floor(start / truncation) + 1) = pos;
            Efield(:, floor(i / truncation) - floor(start / truncation) + 1) = E_grid;
        end
    end

    % plot in real time 
    if plotRealTime || (i == Nt)
        clf; % Clear current figure
        if plotElectricField && plotPos
            subplot(3,1,1);
            histogram(vel, 300);
            xlim([-8 8]);
            xlabel('v');
            ylabel('counts');
            title('Velocity Distribution');
            
            subplot(3,1,2);
            histogram(pos, 300);
            xlim([0 boxsize]);
            xlabel('x');
            ylabel('counts');
            title('Position Distribution');

            subplot(3,1,3);
            plot(linspace(0, boxsize, Nx), E_grid);
            xlim([0 boxsize]);
            xlabel('x');
            ylabel('Electric Field (E)');
            title('Electric Field');
        elseif plotElectricField
            subplot(2,1,1);
            histogram(vel, 300);
            xlim([-8 8]);
            xlabel('v');
            ylabel('counts');
            title('Velocity Distribution');
            
            subplot(2,1,2);
            plot(linspace(0, boxsize, Nx), E_grid);
            xlim([0 boxsize]);
            xlabel('x');
            ylabel('Electric Field (E)');
            title('Electric Field');
        elseif plotPos
            subplot(2,1,1);
            histogram(vel, 300);
            xlim([-8 8]);
            xlabel('v');
            ylabel('counts');
            title('Velocity Distribution');
            
            subplot(2,1,2);
            histogram(pos, 300);
            xlim([0 boxsize]);
            xlabel('x');
            ylabel('counts');
            title('Position Distribution');
        else
            histogram(vel, 300);
            xlim([-8 8]);
            xlabel('v');
            ylabel('counts');
            title('Velocity Distribution');
        end
        drawnow;
    end
end

if storeValues
    % Assign values to global workspace
    assignin('base', 'Vel', Vel);
    assignin('base', 'Pos', Pos);
    assignin('base', 'Efield', Efield);
end


