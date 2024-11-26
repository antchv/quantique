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
w = 1;

f = @(x,n) VA + (VB - VA) * exp(-(2*x/w).^(2*n));
pqrec = @(x) VA*(abs(x)>.5);

parfor p=1:length(NN),p;   % parallel for

  N = floor(NN(p));

  delta = L_bar/N;
  x_bar = -L_bar/2 + L_bar/N * (0:N);

  %vn = pqrec(x_bar);
  vn = f(x_bar,10);

  ee = ones(N+1,1);

  Lap = spdiags([-ee/12 4*ee/3 -5*ee/2 4*ee/3 -ee/12],[-2 -1 0 1 2],N+1,N+1);    %spdiags = sp pour sparse (matrice creuse)
  
  A = -1 /(pi^2 * delta^2) * Lap + spdiags(vn.',0,N+1,N+1);
  options = struct('tol', 1e-6, 'maxit', 1000);
  [psi, En] = eigs(A, n_modes, 'sm', options);

  EEn(:,p) = E0 * sort(diag(En));
  
end


semilogx(NN,EEn,'Linewidth',2);

