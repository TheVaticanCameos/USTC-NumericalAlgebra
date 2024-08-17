a = 0.5;
eps = 1;
n = 100;
h = 1/n;
x = 0:h:1;
yExact = (1-a)*(1-exp(-x/eps))/(1-exp(-1/eps)) + a*x;
