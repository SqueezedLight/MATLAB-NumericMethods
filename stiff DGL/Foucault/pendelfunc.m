function f = pendelfunc(t,w)

L = 67;
g = 9.83;
lambda = 49*pi/180;
Omega = 2*pi/86400;

f = zeros(4,1);
f(1) = w(3);
f(2) = w(4);
f(3) = 2*Omega*sin(lambda)*w(4) - (g/L)*w(1);
f(4) = -2*Omega*sin(lambda)*w(3) - (g/L)*w(2);