clear all;
clc;

t_exact = linspace(0,4,100); 
y_exact1 = 3 - 2*exp(-t_exact);  %exakte ergebnisse
y_exact2 = 3 - 0.998*exp(-1000*t_exact) - 2.002*exp(-t_exact);

time = [0 4];   %Intervallgrenzen
y_01 = 1;        %Anfangsbedingung
y_02 = 0;

TOL = 10^-6;    %Genauigkeit d. num. Berechnung
options = odeset('RelTol',TOL);


%% Test 1, ordinary
[tout,y] = ode45('odefunc1',time,y_01,options);

figure('Name','ode45, ordinary','NumberTitle','off');
plot(t_exact,y_exact1,'-r',tout,y,'ob');
xlabel('t');
ylabel('y');
ylim([0 3]);
title(['benoetigte Stuetzpunkte: ', num2str(numel(tout))]);


%% Test 2, stiff
[tout,y] = ode45('odefunc2',time,y_02,options);

figure('Name','ode45, stiff','NumberTitle','off');
plot(t_exact,y_exact2,'-r',tout,y,'ob');
xlabel('t');
ylabel('y');
title(['benoetigte Stuetzpunkte: ', num2str(numel(tout))]);   