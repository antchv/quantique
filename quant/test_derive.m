clear;
nb_pt = 1e3;
delta = 10/nb_pt; 
xi = -5+(0:nb_pt)*delta; 


n=10;

f=@(x) 1-exp(-x.^n);

f_prime = @(x) n*x.^(n-1).*exp(-x.^(n));

DD = (f(xi+delta)-f(xi))/delta;

DG = (f(xi)-f(xi-delta))/delta;

DC = (f(xi+delta)-f(xi-delta))/delta/2;


DS3 = (f(xi-delta) - 2*f(xi)+f(xi+delta))/delta**2;

DS5 = (-f(xi-2*delta) + 16*f(xi-delta) -30*f(xi) + 16*f(xi+delta) -f(xi+2*delta)) / (delta**2)/12;

subplot(311);plot(xi, DD, xi, DG, xi, DC, xi, f_prime(xi));

subplot(312);plot(xi, DS3, xi, DS5);