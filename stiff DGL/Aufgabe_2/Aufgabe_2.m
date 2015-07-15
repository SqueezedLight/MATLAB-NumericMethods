clear all;
clc;

t_exact = linspace(0,4,100);
time = [0 4];

%% Test 1 des Euler-Programmes mit vorgegebenen h und gewöhnlicher DGL
% y_exact = 3 - 2*exp(-t_exact);    %exakte Loesungsfunktion
% 
% h = [0.8,0.5,0.25,0.1];   %definition der schrittweiten
% 
% for run=1:numel(h)
% 
% i_max = (time(end) - time(1))/h(run);
% y_dach = zeros(i_max,1);
% 
% y_dach(1) = 1;
% t = [0,h(run):h(run):4];
% 
% for i=1:i_max
% 
%     y_dach(i+1) = y_dach(i) + h(run)*(3 - y_dach(i));
% end
% 
% 
% if run == 1, figure('Name','Test 1 des Euler-Programmes mit vorgegebenen h','NumberTitle','off'); end
% subplot(2,2,run);
% plot(t_exact,y_exact,'-r',t,y_dach,'ob');
% title(['h = ',num2str(h(run))]);
% xlabel('t');
% ylabel('y');
% 
% end


%% Test 2 des Euler-Programmes mit vorgegebenen h und steifer DGL
y_exact = 3 - 0.998*exp(-1000*t_exact) - 2.002*exp(-t_exact);

h = [0.0010,0.0012,0.0014,0.0016,0.0018,0.0020];

for run=1:6

i_max = floor((time(end) - time(1))/h(run));
y_dach = zeros(i_max,1);

y_dach(1) = 0;
t = [0,h(run):h(run):4];

for i=1:i_max

    y_dach(i+1) = y_dach(i) + h(run)*(-1000*y_dach(i) + 3000 - 2000*exp(-t(i)));
end


if run == 1, figure('Name','Test 2','NumberTitle','off'); end
subplot(2,3,run);
plot(t_exact,y_exact,'-r',t,y_dach,'ob');
title(['h = ',num2str(h(run))]);
xlabel('t');
ylabel('y');

end


%% Test 2 Idee
% y_exact = 3 - 0.998*exp(-1000*t_exact) - 2.002*exp(-t_exact);
% 
% h1 = 0.0015;
% h2 = 0.0030;
% 
% 
% i_max1 = floor((2 - 0)/h1);
% i_max2 = floor((4 - 2)/h2);
% y_dach1 = zeros(i_max1,1);
% y_dach2 = zeros(i_max2,1);
% 
% y_dach1(1) = 0;
% t1 = [0,h1:h1:2];
% t2 = [2,2+h2:h2:4];
% 
% for i=1:i_max1
%     y_dach1(i+1) = y_dach1(i) + h1*(-1000*y_dach1(i) + 3000 - 2000*exp(-t1(i)));
% end
% 
% 
% %% Plot Test 2 mit Idee, erster Bereich
% figure('Name','Test 2 mit Idee, erster Bereich','NumberTitle','off');
% 
% 
% plot(t_exact,y_exact,'-r',t1,y_dach1,'ob');
% title(['h1 = ',num2str(h1)]);
% xlim([0 2]);
% xlabel('t');
% ylabel('y');
% 
% 
% %% weiter für Test 2 mit Idee
% 
% y_dach2(1) = y_dach1(end);
% 
% for i=1:i_max2
%     y_dach2(i+1) = y_dach2(i) + h2*(-1000*y_dach2(i) + 3000 - 2000*exp(-t2(i)));
% end
% 
% %% Plot Test 2 mit Idee, zweiter Bereich
% figure('Name','Test 2 mit Idee, zweiter Bereich','NumberTitle','off');
% 
% 
% plot(t2,y_dach2,'ob');
% title(['h2 = ',num2str(h2)]);
% xlim([2 4]);
% xlabel('t');
% ylabel('y');
% 
% %% Plot Test 2 mit Idee, gesamt
% figure('Name','Test 2 mit Idee, gesamt','NumberTitle','off');
% 
% 
% plot(t1,y_dach1,'ob',t2,y_dach2,'ob');
% xlabel('t');
% ylabel('y');

