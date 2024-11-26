%------------------- Le main -----------------------
clear;
a = 10e-9;          % largeur d'un puits quantique en m
Npq = 10;           % nombre de puits quantiques
Nmodes = 50;        % nombre de modes à extraire
N = 1e4;            % nombre de points de discrétisation
V0 = 400;           % hauteur de la barrière de potentiel en meV
d = linspace(.1,10,100)*1e-9;  % distance inter-puits en m

% Initialisation d'une matrice pour stocker les niveaux d'énergie
En = zeros(Nmodes, length(d));

% Boucle sur différentes distances inter-puits
for p = 1:length(d)
    En(:, p) = MultiPQs(a, d(p), Npq, V0, Nmodes, N);
end

% Tracer des résultats
figure;
plot(d * 1e9, En, 'Linewidth', 2);
aa = 20; grid; set(gca, 'FontSize', aa);
xlabel('d (nm)', 'FontSize', aa + 10, 'FontWeight', 'normal', 'FontName', 'Times', 'Color', 'b');
ylabel('{ E_n}', 'FontSize', aa + 10, 'FontWeight', 'normal', 'Color', 'b')
%---------------------------------------------