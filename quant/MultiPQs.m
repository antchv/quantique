%-------------- La fonction ----------------
function En = MultiPQs(a, d, Npq, V0, Nmodes, N)
    % Paramètres physiques
    a0 = a;
    me = 9.1091e-31; % masse de l'électron en kg
    meff = 0.067 * me; % masse effective
    e = 1.602176565e-19; % charge élémentaire en C
    hbar = 6.626e-34 / (2 * pi); % constante de Planck réduite en J.s

    % Énergie de seuil pour normalisation
    E0 = hbar^2 * pi^2 / (2 * meff * a0^2) / e * 1e3;

    % Mode 1 pour un puits infini
    V0b = V0 / E0; % normalisation de la hauteur de la barrière
    options.disp = 0;

    % Création de la matrice de potentiel
    ab = a / a0;
    db = d / a0;
    xbm = -10;
    xbM = 10 + Npq * ab + (Npq - 1) * db;
    delt = (xbM - xbm) / N;
    xb = xbm + delt * (0:N);
    f = @(xb, ab, db) (xb >= 0 & xb <= ab);
    
    % Construction du potentiel pour Npq puits quantiques
    for pp = 2:Npq
        f = @(xb, ab, db) f(xb, ab, db) + (xb >= (pp - 1) * (ab + db) & xb <= (pp - 1) * (ab + db) + ab);
    end
    vn = V0b * (1 - f(xb, ab, db));

    % Affichage du potentiel
    plot(xb, vn, 'Linewidth', 2);
    aa = 20; grid;
    set(gca, 'FontSize', aa);
    pause(.01)
    ylabel('{\it \bar{x}}', 'FontSize', aa + 10, 'FontWeight', 'normal', 'Color', 'b')

    % Création de la matrice Hamiltonienne
    ee = ones(N + 1, 1);
    Lap = spdiags([ee -2 * ee ee], [-1 0 1], N + 1, N + 1);
    A = -1/pi^2/delt^2 * Lap + spdiags(vn.', 0, N + 1, N + 1);
    
    % Résolution du problème aux valeurs propres
    [psi, En] = eigs(A, Nmodes, 'sm', options);
    En = E0 * diag(En);
    
    % Suppression des valeurs d'énergie au-dessus de la barrière
    En(En > V0) = nan;
end
%---------------------------------------------


