function En = function_pqrecfull(a, d, V0, N, Lb, modes, p)
    % This function calculates the eigenvalues (energies) of a quantum well
    % potential defined by a series of rectangular potentials.

    % Constants related to the electron
    m_electron = 9.1091e-31;
    e = 1.60217e-19;

    % Effective mass and reduced Planck constant
    m_eff = 0.067 * m_electron;
    hbar = 6.626e-34 / (2 * pi);

    % Reference energy E0 and conversion to meV
    E0 = hbar^2 * pi^2 / (2 * m_eff * a^2);
    E0 = E0 * 1000 / e;

    % Depth of potential well
    v0 = V0 / E0;

    % Discretization of space
    x = -Lb/4 + Lb/N * (0:N);

    % Initialization of the rectangular potential vector
    PQ_rectangulaire = zeros(1, N + 1);



    % Initialization of the triangular potential vector
    % V_triangular = zeros(1, N + 1);

    % % Loop to add triangular potentials
    % for i = 1:num_wells
    %     well_start = -Lb/4 + (i-1) * Lb/(num_wells);
    %     well_end = -Lb/4 + i * Lb/(num_wells);
    %     V_triangular = V_triangular + epsilon * e * x .* (x >= well_start & x <= well_end);
    % end

    % % Normalizing the potential
    % V_triangular = V_triangular - min(V_triangular);



    % Loop to add rectangular potentials
    for i = 1:p
        PQ_rectangulaire = PQ_rectangulaire + (x >= (i - 1) * (1 + d/a) & x <= (i + (i - 1) * d/a));
    end

    % Normalizing the potential
    PQ_rectangulaire = v0 - v0 * PQ_rectangulaire;

    % Discretization step
    delt = Lb / N;

    % Potential vector
    vn = PQ_rectangulaire;

    % Vector of ones
    ee = ones(N + 1, 1);

    % Laplacian matrix using finite differences
    Lap = spdiags([-ee/12 4*ee/3 -5*ee/2 4*ee/3 -ee/12], [-2 -1 0 1 2], N + 1, N + 1);

    % Diagonal matrix of the potential
    vi = spdiags(vn.', 0, N + 1, N + 1);

    % Hamiltonian matrix A of the Schrödinger equation
    A = -1/(pi*delt)^2 * Lap + vi;

    % Calculating eigenvalues and eigenvectors using the 'eigs' function
    %%%% RETOURNER PSI POUR AVOIR LES FONCTIONS D'ONDE ASSOCIE AU ENERGIE

    [psi, En] = eigs(A, modes, 'sm');

    % Sorting and converting eigenvalues to meV
    En = E0 * sort(diag(En));

    % Plotting the potential profile
    subplot(2, 1, 2);
    plot(x, PQ_rectangulaire, 'LineWidth', 2);




   
    % Plotting the potential profile
    % subplot(2, 1, 2);
    % plot(x, PQ_rectangulaire, 'LineWidth', 2);
    % hold on;

    % % Plotting the wavefunctions
    % for j = 1:modes
    %     plot(x, psi(:, j) + En(j), 'LineWidth', 1.5);
    % end
    % hold off;
    
end
