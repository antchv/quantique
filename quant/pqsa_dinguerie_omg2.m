% COnvergence cas semi-analytique

% clear


% A = linspace(1e-9,10e-9,500);    % ordre de grandeur de la taille du puit quantique (ici 10 nm)
% V0 = 1e3;
% N=1000;
% n_modes = 3;
% L_bar = 5;
% EEn1=zeros(length(A));
% EEn2=zeros(length(A));
% EEn3=zeros(length(A));
% parfor a=1:length(A)
%   E=function_pqrec(A(a), V0, N, L_bar, n_modes);
%   EEn1(a)=E(1);
%   EEn2(a)=E(2);
%   EEn3(a)=E(3);
% end

% plot(A, EEn1, A, EEn2, A, EEn3);

% Convergence cas semi-analytique

clear

A = linspace(1e-10, 10e-9, 100);  % ordre de grandeur de la taille du puit quantique (ici 10 nm)
V0 = 400;
N = 1000;
n_modes = 4;
L_bar = 5;

% Préallocation de la matrice des résultats
EEn = zeros(length(A), n_modes);

% Précalcul des valeurs constantes si possible

parfor a = 1:length(A)
    E = function_pqrec(A(a), V0, N, L_bar, n_modes);
    EEn(a, :) = E;
end

% Plot pour chaque mode
figure;
hold on;
for mode = 1:n_modes
    plot(A*1e9, EEn(:, mode), 'DisplayName', sprintf('Mode %d', mode));
end
hold off;
legend;
grid();
axis([0 10 0 400])

 

