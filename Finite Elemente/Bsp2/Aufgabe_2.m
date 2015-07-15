clear all;
clc;

A = textread('./dreiecke_126.txt');

n = A(1,1);
npkt = A(1,2);
cborder = [1,5,6,7,8,9,10,11,2];
dborder = [4,27,26,25,24,23,22,21,3];
wert = [30,30,40,60,75,60,50,30,50];

T_inter = zeros(101,101);
c = [1,0,0;-1,1,0;-1,0,1];

a0 = 2;
a1 = 2;
a2 = 2;
a4 = [-0.3,-0.3,0.2,-0.3];
a5 = 0;
dt = 0.5;
tschritte = 100;

x = A(3+n:end,1);
y = A(3+n:end,2);
xnetz = linspace(0,4,101);
ynetz = linspace(0,4,101);
M2 = [2,1,1;1,2,1;1,1,2];

for viermal=1:4
inhom = zeros(npkt,1);
mata = zeros(npkt,npkt,n);
matb = zeros(npkt,npkt,n);
B = zeros(npkt,npkt,n);
h = [0,10,30,-10];
stat = [103.626,141.290,75.244,65.961];

for m=1:n,    
    i = A(2+m,1);
    j = A(2+m,2);
    k = A(2+m,3);
    ijk = [i,j,k];
    
    d = (x(j)-x(i))*(y(k)-y(i)) - (x(k)-x(i))*(y(j)-y(i));      %Jakobi Determinante
    
    e11 = (a1*(y(k)-y(i))^2 + a2*(x(k)-x(i))^2)/d;
    e12 = -(a1*(y(k)-y(i))*(y(j)-y(i)) + a2*(x(k)-x(i))*(x(j)-x(i)))/d;
    e22 = (a1*(y(j)-y(i))^2 + a2*(x(j)-x(i))^2)/d;
    
    M1 = [e11+2*e12+e22,-e11-e12,-e12-e22;-e11-e12,e11,e12;-e12-e22,e12,e22];
    b = h(viermal)*d/6.*[1;1;1];
    ma = 1/2*M1;
    Bmat = (a0*d/24).*M2;
    
    for lauf=1:3,
        for lauff=1:3
            %% Beitrag der einzelnen Dreiecke
            B(ijk(lauff),ijk(lauf),m) = Bmat(lauff,lauf);
            mata(ijk(lauff),ijk(lauf),m) = ma(lauff,lauf);
            inhom(ijk(lauff),m) = b(lauf);
    
            %% Beitrag aus den Causchy'schen Randbedingungen
            if isempty(cborder(cborder==i)) == 0 && isempty(cborder(cborder==j)) == 0 && lauf<3 && lauff<3,
                d_pq = sqrt((x(i)-x(j))^2 + (y(i)-y(j))^2);
                mb = a4(viermal)*d_pq/6.*[2,1;1,2];
                matb(ijk(lauff),ijk(lauf),m) = mb(lauff,lauf);
            end
        end
    end
    
end

%% Über alle Durchgänge summieren
mat = sum(mata + matb,3);  %mata...Beitrag Dreiecke   matb...Beitrag Cauchy
inhom = sum(inhom,2);
B = sum(B,3);
C = mat + (2/dt).*B;
D = mat - (2/dt).*B;

%% T_anf berechnen
T_anf = 100.*ones(npkt,1);
T_anf(dborder) = wert;

T_old = T_anf;
    
%% Dirichlet berücksichtigen
iwert = 1;

for ind=dborder,  
    for np=1:npkt,
       inhom(np) = inhom(np) - mat(np,ind)*wert(iwert);
       C(ind,np) = 0;
       C(np,ind) = 0;
       D(ind,np) = 0;
       D(np,ind) = 0;
    end
    iwert = iwert + 1;
end

iwert = 1;

for ind=dborder, 
    C(ind,ind) = 1;
    D(ind,ind) = -1;
    inhom(ind) = 0;
    T_old(ind) = wert(iwert);
    iwert = iwert + 1;
end

    
%% Schleife über die Zeit
for z=1:tschritte+1,
    T_dt = C\(2.*inhom - D*T_old);
    T_old = T_dt;
    
    for lauff=1:101
        for lauf=1:101,      
                for m=1:n,    
                    i = A(2+m,1);
                    j = A(2+m,2);
                    k = A(2+m,3);
                    
                    v = (xnetz(lauf)*(y(i) - y(j)) + ynetz(lauff)*(x(j) - x(i)) - y(i)*x(j) + x(i)*y(j))/(y(k)*(x(j)-x(i)) + y(i)*(x(k) - x(j)) + y(j)*(x(i) - x(k)));
                    u = (xnetz(lauf) - x(i) - (x(k) - x(i))*v)/(x(j) - x(i));

                    if v<=-u+1 && u>=0 && u<=1 && v>=0,
                        cm = c*[T_dt(i); T_dt(j); T_dt(k)];
                        T_inter(lauff,lauf) = cm(1) + cm(2)*u + cm(3)*v;
                    end
                end
        end    
    end
    
%     subplot(2,2,viermal);
%     pcolor(T_inter);
%     shading flat;
%     caxis([0,150]);
%     axis('equal');
%     xlabel('x-Achse');
%     ylabel('y-Achse');
%     title(['Animation: t = ',num2str(z*dt)]);
%     M(z) = getframe;

%  if viermal == 1 && z == 1, q = figure; 
%     else figure(q);
%     end
    hold on;
    subplot(2,2,viermal);
    pcolor(T_inter);
    shading flat;
    caxis([0,150]);
    axis('equal');
    xlabel('x-Achse');
    ylabel('y-Achse');
    title(['Animation: t = ',num2str((z-1)*dt),'     Durchlauf #',num2str(viermal)]);
    M(z) = getframe;
    hold off;

%     if viermal == 1 && z == 1, g = figure; 
%     else figure(g);
%     end
%     hold on;
%     subplot(2,2,viermal);
%     plot(z,T_dt(2),z,stat(viermal));
%     xlabel('Zeitschritte [s]');
%     ylabel('Temperatur [°C]');
%     title(['Durchlauf #', num2str(viermal)]);
%     %legend('Temperaturverlauf T2', 'Stationäre Lösung');
%     xlim([0,100]);
%     hold off;
end

%movie(M,1);

end