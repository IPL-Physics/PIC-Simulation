%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the simulation
%Na,Nb, and Nc are the counts of three separate distributions of
%electrons. tEnd is how long the simulation will go for. va,vb, and vc
%are beam velocities or mean velocity of each of the separate
%distributions. and vtha,vthb, and vthc are the variance or thermal
%velocities of those distributions
%plasma_pic_simulation(particleA,particleB,particleC,tEnd,vA,vB,vC,TA,TB,TC);
plasma_pic_simulation(30000,9000,0,200,0,4,0,0.5,0.2,0);

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Phase Space Plots %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

timestep = size(Vel, 2);
num_bins = [300, 300]; % Adjust the number of bins as needed

% Create figure
figure;
h = histogram2(Pos(:,1), Vel(:,1), 'DisplayStyle', 'tile','ShowEmptyBins', 'off', 'NumBins', num_bins);
colorbar;
xlabel('Position');
ylabel('Velocity');
zlabel('Counts');
title('Phase Space Plot');
colorbar; % Add colorbar
view(2); % Set view to 2D

% Animation loop
for t = 1:timestep
    % Calculate histogram data
    [N, edgesX, edgesY] = histcounts2(Pos(:,t), Vel(:,t), num_bins);
    
    % Delete the old histogram if it exists
    if exist('h', 'var')
        delete(h);
    end
    
    % Create a new histogram2 plot
    h = histogram2(Pos(:,t), Vel(:,t), num_bins);
    xlabel('Position');
    ylabel('Velocity');
    zlabel('Frequency');
    title(sprintf('Phase Space Plot - Time Step %d', t));
    colorbar; % Add colorbar
    view(2); % Set view to 2D

    % Set the bin edges to match the new histogram
    h.XBinEdges = edgesX;
    h.YBinEdges = edgesY;
    
    % Draw updated plot
    drawnow;
    
    % Pause for a short period to visualize the animation (adjust as needed)
    pause(0.0);
end

%% 
% Calculate mean, variance, skewness, and kurtosis for each column
means = mean(Vel, 1);
variances = var(Vel, 0, 1); % The second argument is set to 0 to get the sample variance
skewness_vals = skewness(Vel, 0, 1);
kurt_vals = kurtosis(Vel, 0, 1);
timestep = size(Vel, 2);

% Create a figure with subplots
figure;
subplot(5, 1, 1);
scatter(1:timestep, means, 'filled', 'MarkerFaceColor', 'blue');
ylabel('Mean');
legend('Mean');

subplot(5, 1, 2);
scatter(1:timestep, variances, 'filled', 'MarkerFaceColor', 'magenta');
ylabel('Variance');
legend('Variance');

subplot(5, 1, 3);
scatter(1:timestep, skewness_vals, 'filled', 'MarkerFaceColor', 'green');
ylabel('Skewness');
legend('Skewness');

subplot(5, 1, 4);
scatter(1:timestep, kurt_vals, 'filled', 'MarkerFaceColor', 'red');
ylabel('Kurtosis');
legend('Kurtosis');

q_array = (6 + 7 * kurt_vals) ./ (6 + 5 * kurt_vals);
subplot(5, 1, 5);
scatter(1:timestep, q_array, 'filled', 'MarkerFaceColor', 'black');
ylabel('nonextensive q');
xlabel('Timestep');
legend('nonextensive q');

% Adjust layout
sgtitle('Mean, Variance, Skewness, Kurtosis, and nonextensive q over Timesteps');
%% 
% Key array sizes
num_timesteps = size(Vel, 2);
num_points = size(Vel, 1);

% Define different smoothing parameters (sigma) for each data type
sigma_Vel = 0.1;
sigma_Pos = 0.2;
sigma_Efield = 0.01;

% Number of bins
num_bins = 300;

% Initialize variables to store smoothed data
Vel_smooth = zeros(num_bins, num_timesteps);
Pos_smooth = zeros(num_bins, num_timesteps);
Efield_smooth = zeros(num_bins, num_timesteps);

% Create a figure and subplots outside the loop
fig = figure;
set(fig, 'Position', get(0, 'ScreenSize')); % Set figure to fullscreen
subplot1 = subplot(4,1,1); % Smoothed Velocities
subplot2 = subplot(4,1,2); % Derivative of Smoothed Velocities
subplot3 = subplot(4,1,3); % Smoothed Positions
subplot4 = subplot(4,1,4); % Smoothed Electric Field

% Loop through each timestep
for t = 1:num_timesteps
    % Create histograms for Vel, Pos, and Efield
    [hist_Vel, bin_edges_Vel] = histcounts(Vel(:, t), num_bins);
    [hist_Pos, bin_edges_Pos] = histcounts(Pos(:, t), num_bins);
    [hist_Efield, bin_edges_Efield] = histcounts(Efield(:, t), num_bins);
    
    % Calculate bin centers
    bin_centers_Vel = (bin_edges_Vel(1:end-1) + bin_edges_Vel(2:end)) / 2;
    bin_centers_Pos = (bin_edges_Pos(1:end-1) + bin_edges_Pos(2:end)) / 2;
    bin_centers_Efield = (bin_edges_Efield(1:end-1) + bin_edges_Efield(2:end)) / 2;
    
    % Define the Gaussian kernels
    dx_Vel = bin_edges_Vel(2) - bin_edges_Vel(1);
    dx_Pos = bin_edges_Pos(2) - bin_edges_Pos(1);
    dx_Efield = bin_edges_Efield(2) - bin_edges_Efield(1);
    
    x_Vel = -3 * sigma_Vel:dx_Vel:3 * sigma_Vel;
    x_Pos = -3 * sigma_Pos:dx_Pos:3 * sigma_Pos;
    x_Efield = -3 * sigma_Efield:dx_Efield:3 * sigma_Efield;
    
    gaussian_Vel = exp(-0.5 * (x_Vel / sigma_Vel).^2);
    gaussian_Pos = exp(-0.5 * (x_Pos / sigma_Pos).^2);
    gaussian_Efield = exp(-0.5 * (x_Efield / sigma_Efield).^2);
    
    % Normalize the Gaussian kernels
    gaussian_Vel = gaussian_Vel / sum(gaussian_Vel);
    gaussian_Pos = gaussian_Pos / sum(gaussian_Pos);
    gaussian_Efield = gaussian_Efield / sum(gaussian_Efield);
    
    % Smooth the histograms using convolution
    smoothed_Vel = conv(hist_Vel, gaussian_Vel, 'same');
    smoothed_Pos = conv(hist_Pos, gaussian_Pos, 'same');
    smoothed_Efield = conv(hist_Efield, gaussian_Efield, 'same');
    
    % Store the smoothed histograms
    Vel_smooth(:, t) = smoothed_Vel;
    Pos_smooth(:, t) = smoothed_Pos;
    Efield_smooth(:, t) = smoothed_Efield;
    
    % Calculate the derivative of the smoothed velocities
    derivative_Vel = diff(smoothed_Vel) / dx_Vel;
    bin_centers_Vel_derivative = bin_centers_Vel(1:end-1) + dx_Vel/2; % Adjust bin centers for derivative
    
    % Update the subplots with the new data
    axes(subplot1);
    plot(bin_centers_Vel, smoothed_Vel);
    title('Smoothed Velocities');
    
    axes(subplot2);
    plot(bin_centers_Vel_derivative, derivative_Vel);
    title('Derivative of Smoothed Velocities');
    xlim([0, 8]);
    ylim([min(derivative_Vel(:)), max(derivative_Vel(:))]);
    % Add horizontal and vertical red lines to the derivative plot
    hold on;
    yline(0, 'r', 'LineWidth', 1); % Horizontal red line at y = 0
    xline(0, 'r', 'LineWidth', 1); % Vertical red line at x = 0
    hold off;

    axes(subplot3);
    plot(bin_centers_Pos, smoothed_Pos);
    title('Smoothed Positions');
    ylim([min(Pos_smooth(:)), max(Pos_smooth(:))]);
    
    axes(subplot4);
    plot(bin_centers_Efield, smoothed_Efield);
    title('Smoothed Electric Field');
    ylim([min(Efield_smooth(:)), max(Efield_smooth(:))]);
    
    % Adjust the figure
    drawnow;
    pause(0.25);
end

%% 

% Compute the FFT and handle sizes
Pos_freq = [];
Efield_freq = [];

num_points = size(Vel, 1);

for t = 1:num_timesteps
    % Compute FFT
    fft_pos = fft(Pos_smooth(:, t));
    fft_pos = fft_pos(1:size(Pos_smooth,1)/2+1); % Only positive frequencies
    fft_efield = fft(Efield_smooth(:, t));
    fft_efield = fft_efield(1:numel(Efield(:,1))/2+1); % Only positive frequencies
    
    % Store the result
    Pos_freq(:, t) = fft_pos;
    Efield_freq(:, t) = fft_efield;
end

% Frequency axis
Fs = 10; % Sampling frequency, adjust as needed
f = (0:(num_points-1))*(Fs/num_points);
fe = (0:(numel(Efield(:,1))-1))*(Fs/numel(Efield(:,1)));

% Create figure for derivative and FFT
figure;
for t = 1:num_timesteps
    
    subplot(2, 1, 1);
    plot(abs(Pos_freq(:, t)),10);
    title('FFT of Smoothed Position');
    xlabel('Frequency');
    ylabel('Magnitude');
    
    subplot(2, 1, 2);
    plot(abs(Efield_freq(:, t)),1000);
    title('FFT of Smoothed Electric Field');
    xlabel('Frequency');
    ylabel('Magnitude');
    
    drawnow;
end

