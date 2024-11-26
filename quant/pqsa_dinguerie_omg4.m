clear

% Parameters
a = 10e-9;         % Width of the potential well
d = linspace(0.1, 5, 100) * 1e-9;  % Array of potential parameters
p = 3;           % Number of potential wells
V0 = 400;         % Depth of potential well
N = 5e3;          % Number of points for spatial discretization
modes = 10;       % Number of eigenmodes to calculate

% Initialize En matrix to store eigenvalues
En = zeros(modes, length(d));

% Loop through different potential parameters
fprintf('Progress: 0%%');
parfor i = 1:length(d)
    % Calculate the effective length Lb
    Lb = 2 * p * (1 + d(i) / a);

    % Call the function to calculate eigenvalues for the current potential
    En(:, i) = function_pqrecfull(a, d(i), V0, N, Lb, modes, p);

    % Call the function to calculate eigenvalues and eigenvectors for display wave function too
    %[En(:, i), ~] = function_pqrecfull(a, d(i), V0, N, Lb, modes, p);

    % Display progress
    progress = i / length(d);
    fprintf('\r\033[1;34mProgress: %.1f%%\033[0m', progress * 100);
end
fprintf('\n'); % Move to the next line after completion

% Plotting the results
xlabel('a (nm)');
ylabel('E_n (meV)');
title('Eigenvalues vs. Width of the Potential Well');
subplot(2, 1, 1);
plot(1e9 * d, En, 'LineWidth', 2);
grid on;
