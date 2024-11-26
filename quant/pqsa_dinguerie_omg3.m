clear

D = linspace(1e-10, 20e-9, 1000);  % ordre de grandeur de la taille du puit quantique (ici 10 nm)
V0 = 400;
a = 6e-9;
N = 1000;
n_modes = 4;
L_bar = 15;

% Préallocation de la matrice des résultats
EEn = zeros(length(D), n_modes);

% Précalcul des valeurs constantes si possible

parfor d = 1:length(D)
    E = function_pqrec2(D(d), a, V0, N, L_bar, n_modes);
    EEn(d, :) = E;
end

% Plot pour chaque mode
figure;
hold on;
for mode = 1:n_modes
    plot(D*1e9, EEn(:, mode), 'DisplayName', sprintf('Mode %d', mode));
end
hold off;
legend;
grid();
axis([0 20 0 400])


 

