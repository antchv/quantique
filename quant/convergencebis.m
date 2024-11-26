clear

% constantes
a = 10e-9;
m_e = 9.1091e-31;
e = 1.60217e-19;
h_bar = 6.626e-34 / 2 / pi;
m_eff = 0.067 * m_e;
E_f = h_bar^2 * pi^2 / (2 * m_eff * a^2); % mode fondamental du PQ infini

E_f = E_f * 1000/e; % conversion en meV



Lb=5;N=[100:100:1e3 2000:1e3:9000];
V0=1e7;v0=V0/E_f;

modes=5;
EEn=[];
for p=1:length(N)
  delt=Lb/N(p);xb=-Lb/2+Lb/N(p)*(0:N(p));vn=v0*(abs(xb)>.5);

  ee=ones(N(p),1);Lap=spdiags([ee -2*ee ee],[-1 0 1], N(p)+1, N(p)+1);vvi=spdiags(vn.',0,N(p)+1, N(p)+1);
  A=-1/pi^2/delt^2*Lap+vvi;
  [psi,En]=eigs(A,modes,'sm');En=diag(En);EEn=[EEn En];
end
plot(N, EEn,'Linewidth',2);
