clear all;
clc;

t_exact = linspace(0,0.3,200);
time = [0 0.3];

%% Gleichungssystem mit impliziter Lösungsformel
y_exact1 = 2*exp(-t_exact) - exp(-1000*t_exact);
y_exact2 = -exp(-t_exact) + exp(-1000*t_exact);

h = [0.0005,0.0017,0.0018,0.0019,0.0020,0.0021];
C = [-998,-1998;999,1999];
y_anf = [1;0];

for run=1:numel(h)

i_max = floor((time(end) - time(1))/h(run));

t = [0,h(run):h(run):0.3];

y_dach = zeros(2,i_max);
y_dach(:,1) = y_anf;

for i=1:i_max
    y_dach(:,i+1) = inv(eye(2) + h(run)*C) * y_dach(:,i);
end

%% Plot der Funktionen (zwischen 1 und 2 waehlen)
if run == 1, figure('Name','Gleichungssystem mit impliziter Lösungsformel, y2','NumberTitle','off'); end
subplot(2,3,run);
hold on;
plot(t_exact,y_exact2,'-r');
plot(t,y_dach(2,:),'ob');
hold off;
title(['Schrittweite: ', num2str(h(run))]);
xlabel('t');
ylabel('y');
xlim([0 0.3]);
ylim([-2 2.5]);

end

