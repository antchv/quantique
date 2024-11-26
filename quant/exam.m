%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%    Antoine CHARVIN    %%
        %% M2 PHYSIQUE NUMERIQUE %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




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

V_0 = 1;
V_0 = 1000 * V_0;
V_0b = V_0 / E_f;

alpha_0 = pi * sqrt(V_0b);

f = @(x) abs(cos(x/2)).*(tan(x/2)>0);


ff = @(x) f(x) - x / alpha_0;


x = linspace(0, 1.5*alpha_0, 1e3);


%%%%%%%%%%%%%%%%%%%%%%%%%%
Lb=5;N=1000;delt=Lb/N;xb=Lb/2+Lb/N*(0:N);vi=V_0b*(abs(xb)>.5);
A=-1/(delt.^2*pi^2);

[psi, Eb] = eig(A);

%% modes paires
disp("energies des modes paires : ")
seed = 2;
mod1 = fsolve(ff, seed);
Eb_mod1 = mod1^2 / pi^2; % 0.75 meV
E_mod1_pair = E_f * Eb_mod1 %42.298 meV

seed = 8;
mod2 = fsolve(ff, seed);
Eb_mod2 = mod2^2 / pi^2; % 6.66 meV
E_mod2_pair = E_f * Eb_mod2 %373.92 meV

seed = 13;
mod3 = fsolve(ff, seed);
Eb_mod3 = mod3^2 / pi^2; 
E_mod3_pair = E_f * Eb_mod3



L=5;
nx = 1000;
xb = linspace(-L/2, L/2, nx);


mod = [mod1, mod2, mod3];
psi_val = zeros(nx, length(mod));

for i=1:length(mod)
  alp = mod(i);
  bet = sqrt(pi^2*V_0b-alp^2);

  Ac=1;
  Bd=2*Ac*exp(bet/2)*cos(alp/2);
  Ag=Bd;

  psi = @(xb) (xb<=-0.5).*(Ag * exp(bet*xb)) + (-0.5<xb & xb<=0.5).*(2*Ac*cos(alp*xb)) + (xb>0.5).*(Bd*exp(-bet*xb));
  psi_val(:,i) = psi(xb);
end
ax(1) = subplot(211);
plot(xb, psi_val);
ax(1) = title("pair");



g = @(x) abs(sin(x/2)).*(cot(x/2)<0);


gg = @(x) g(x) - x / alpha_0;


x = linspace(0, 1.5*alpha_0, 1e3);


%%% modes impaires
disp("energies des modes impaires : ")
seed = 5;
mod1 = fsolve(gg, seed);
Eb_mod1 = mod1^2 / pi^2; 
E_mod1_impair = E_f * Eb_mod1

seed = 12;
mod2 = fsolve(gg, seed);
Eb_mod2 = mod2^2 / pi^2;
E_mod2_impair = E_f * Eb_mod2





L=5;
nx = 1000;
xb = linspace(-L/2, L/2, nx);

mod = [mod1, mod2];
psi_val = zeros(nx, length(mod));

for i=1:length(mod)
  alp = mod(i);
  bet = sqrt(pi^2*V_0b-alp^2);

  Ac=1j;
  Bd=2j*Ac*exp(bet/2)*sin(alp/2);
  Ag=-Bd;

  psi = @(xb) (xb<=-0.5).*(Ag * exp(bet*xb)) + (-0.5<xb & xb<=0.5).*(2j*Ac*sin(alp*xb)) + (xb>0.5).*(Bd*exp(-bet*xb));
  psi_val(:,i) = psi(xb);
end
ax(2)=subplot(212);

plot(xb, psi_val);
ax(2)=title("impair");