clear;

va = 1;
vb = 0;
w = 2;

x = linspace(-2*w, 2*w, 1e3);
n = linspace(1, 10, 10);

f = @(x,n) va + (vb-va) * exp(-(2*x/w).^(2*n));

for i = 1:length(n)
  plot(x, f(x,n(i)));
  hold on;
end
hold off
