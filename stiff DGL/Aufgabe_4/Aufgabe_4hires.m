clear all;
clc;

time = [0 400];   %Intervallgrenzen

y0 = zeros(8,1);
y0(1) = 1;        %Anfangsbedingung
y0(8) = 0.0057;

TOL = 5*10^-4;    %Genauigkeit d. num. Berechnung


%% HIRES

options = odeset('RelTol',TOL,'InitialStep',0.001);

[tout,y] = ode45('aufgabe4hiresfunc',time,y0,options);

h = figure('Name','HIRES Problem mit ode23s','NumberTitle','off');


for i=1:8
    subplot(2,4,i);
    plot(tout,y(:,i));
    
    if i<5 || i>6, xlim([0 5]);
    else xlim([0 400]); end
    
    xlabel('t');
    ylabel(['y',num2str(i)]);
end

       
        
        