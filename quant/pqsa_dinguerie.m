clear
y=[0,1];
% constantes
a = 10e-9;
m_e = 9.1091e-31;
e = 1.60217e-19;
h_bar = 6.626e-34 / 2 / pi;
m_eff = 0.067 * m_e;
E_f = h_bar^2 * pi^2 / (2 * m_eff * a^2); % mode fondamental du PQ infini

E_f = E_f * 1000/e; % conversion en meV


V_0 = 1e15 ;
Vb = V_0 / E_f;

Lb=5;N=10000;delt=Lb/N;xb=-Lb/2+Lb/N*(0:N);vi=Vb*(abs(xb)>.5);

lap = -2*diag(ones(1,N)) + diag(ones(1, N-1), 1) + diag(ones(1, N-1), -1);
modes=3;
ee=ones(N+1,1);Lap=spdiags([ee -2*ee ee],[-1 0 1], N+1, N+1);vvi=spdiags(vi.',0,N+1,N+1);


A=-Lap./(delt.^2*pi^2) + vvi;


[psi, Eb] = eigs(A, modes,'sm');

#Eb=E_f*diag(Eb);
Eb=diag(Eb);

Eb;


EEn=[];
V0=logspace(3,7,100);v0=V0/E_f;

parfor p=1:length(V0)
  vn=v0(p)*(abs(xb)>.5);
  A=1/pi^2/delt^2*Lap+spdiags(vn.',0,N+1,N+1);
  [psi,En]=eigs(A,modes,'sm');En=diag(En);EEn=[EEn En];
end
#plot(xb, psi(:,1:1:2));
#
#subplot(321);plot(xb, psi(:,7));
#subplot(322);plot(xb, psi(:,6));
#subplot(323);plot(xb, psi(:,5));
#subplot(324);plot(xb, psi(:,4));
#subplot(325);plot(xb, psi(:,3));
#subplot(326);plot(xb, psi(:,2));
plot(v0, EEn);
