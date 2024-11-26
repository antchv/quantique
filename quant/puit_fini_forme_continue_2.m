% COnvergence cas semi-analytique

clear

a = 10e-9;    % ordre de grandeur de la taille du puit quantique (ici 10 nm)
m_e = 9.1091e-31;
m_eff = .067*m_e;
e = 1.602176565e-19;
h_bar = 6.626e-34/(2*pi);   % m_eff dans l' AsGa (1.08

E0 = h_bar^2*pi^2/(2*m_eff*a^2)/e*1e3;  % Mode 1 pour le puit infini en meV
%V0 = [1e3 1e4 1e5 1e6 1e7 1e8];
V0 = 1e3;
v0 = V0/E0;   % V0 en meV puis normalisation

n_modes = 2;
options.disp = 0;

L_bar = 5;
NN = logspace(2,5,50);

VA = v0;
VB = 0;
w = 1/2;

f = @(x,n) VA + (VB - VA) * exp(-(2*x/w).^(2*n));

parfor p=1:length(NN),p   % parallel for

  N = floor(NN(p));

  delta = L_bar/N;
  x_bar = -L_bar/2 + L_bar/N * (0:N);

  vn = v0 * (abs(x_bar)>.5);
  %vn = f(x_bar,4)

  ee = ones(N+1,1);
  Lap = spdiags([ee -2*ee ee],[-1 0 1],N+1,N+1);    %spdiags = sp pour sparse (matrice creuse)

  A = -1 /(pi^2 * delta^2) * Lap + spdiags(vn.',0,N+1,N+1);
  [psi,En] = eigs(A,n_modes,'sm',options);

  EEn(:,p) = E0 * sort(diag(En));
end

% Energie du mode n cas Semi-Analytique
E1SA = 42.298056780211112; %meV
E2SA = 168.1689175719088; %meV

E_vrai = [E1SA;E2SA];
E_vrai = repmat(E_vrai,1,length(NN));

semilogx(NN,abs(EEn-E_vrai),'Linewidth',2)
