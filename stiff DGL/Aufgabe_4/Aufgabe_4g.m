clear all;
clc;

time = [0 0.5];   %Intervallgrenzen
t_exact = linspace(0,0.5);

y1_exact = 2*exp(-t_exact) - exp(-1000*t_exact);
y2_exact = -exp(-t_exact) + exp(-1000*t_exact);

y0 = zeros(2,1);
y0(1) = 1;        %Anfangsbedingung

TOL = 10^-4;    %Genauigkeit d. num. Berechnung


%% Gear Problem

options = odeset('RelTol',TOL);

[tout,y] = ode23s('aufgabe4gfunc',time,y0,options);

h = figure('Name','Gear Problem mit ode23s','NumberTitle','off');
plot(t_exact,y1_exact,'-r',tout,y(:,1),'ob',t_exact,y2_exact,'-g',tout,y(:,2),'ob');
legend('exakt y1','num. y1','exakt y2','num y2');
title(['TOL = ', num2str(TOL)]);
xlabel('t');
ylabel('y');
 

       
        
        