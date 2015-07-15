clear all;
clc;

L = 67;
g = 9.83;
lambda = 49*pi/180;
Omega = 2*pi/86400;

time = [0 43200];   %Intervallgrenzen
w0 = zeros(4,1);
w0(1) = L/100;        %Anfangsbedingung

TOL = [10^-2,10^-3,10^-4,10^-5,10^-6];    %Genauigkeit d. num. Berechnung
enter1 = 1;
enter2 = 1;
E0 = g*(L-sqrt(L^2-(L/100)^2));

dat1 = zeros(4,3);
dat2 = zeros(3,3);

%% Foucault'sches Pendel

for i=1:5,
    options = odeset('RelTol',TOL(i));
    
    
% %% ode23 (allgemein anwendbare Routine, explizites Rechenverfahren)
%     if 2 <= i && i <= 5,
%         tic;
%         [tout,w] = ode23('pendelfunc',time,w0,options);
%         tElapsed=toc;
%         if enter1 == 1, b = figure('Name','Foucaultsches Pendel mit explizitem Verfahren (ode23)','NumberTitle','off'); end
%         figure(b);
%         subplot(2,2,enter1);
%         axis square;
%         plot(w(:,1),w(:,2));
%         title(['TOL = ',num2str(TOL(i))]);
%         xlabel('x');
%         ylabel('y');
%         
%         %% Berechnung der Gesamtenergie
%         Et = (w(:,3).^2 + w(:,4).^2)/2 + g.*(L - sqrt(L^2 - w(:,1).^2 - w(:,2).^2));
%         if enter1 == 1, c = figure('Name','Normierte Gesamtenergie, ode23','NumberTitle','off'); end
%         figure(c);
%         subplot(2,2,enter1);
%         plot(tout,Et / E0);
%         title(['TOL = ',num2str(TOL(i))]);
%         xlabel('t');
%         ylabel('Et/E0');
%         
%         enter1 = enter1 + 1;
% 
%         %% Tabelle erstellen
%         dat1(i-1,:) = [TOL(i),numel(tout),tElapsed];
%         
%         if i == 2, l = figure('Position',[100 100 350 150],'Name','Tabelle ode23','NumberTitle','off'); end
%         figure(l);
%         cnames = {'TOL','Stuezpkt.','Rechenzeit'};
%         t = uitable('Data',dat1,'ColumnName',cnames,'Parent',l,'Position',[20 20 300 100]);
%         
%     end

%% ode23t (auf stiff DGL spezialisiert, implizit)
    if 1 <= i && i <= 3,
        tic;
        [tout,w] = ode23t('pendelfunc',time,w0,options);
        tElapsed=toc;
        if enter2 == 1, d = figure('Name','Foucaultsches Pendel mit implizitem Verfahren (ode23t)','NumberTitle','off'); end
        figure(d);
        subplot(2,2,enter2);
        axis square;
        plot(w(:,1),w(:,2));
        title(['TOL = ',num2str(TOL(i))]);
        xlabel('x');
        ylabel('y');
        
        %% Berechnung der Gesamtenergie
        Et = (w(:,3).^2 + w(:,4).^2)/2 + g.*(L - sqrt(L^2 - w(:,1).^2 - w(:,2).^2));
        if enter2 == 1, e = figure('Name','Normierte Gesamtenergie, ode23t','NumberTitle','off'); end
        figure(e);
        subplot(2,2,enter2);
        plot(tout,Et / E0);
        title(['TOL = ',num2str(TOL(i))]);
        xlabel('t');
        ylabel('Et/E0');

        enter2 = enter2 + 1;
        
        %% Tabelle erstellen
        dat2(i,:) = [TOL(i),numel(tout),tElapsed];
        
        if i == 1, k = figure('Position',[100 100 350 150],'Name','Tabelle ode23t','NumberTitle','off'); end
        figure(k);
        cnames = {'TOL','Stuezpkt.','Rechenzeit'};
        t = uitable('Data',dat2,'ColumnName',cnames,'Parent',k,'Position',[20 20 300 100]);
        

    end
    


end




