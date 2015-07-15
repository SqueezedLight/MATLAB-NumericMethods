clear all;
clc;

A = textread('./Vieleck_680.txt');

n = A(1,1);
npkt = A(1,2);

h = 0;
a0 = 0;
a1 = 1;
a2 = 1;
a4 = 0;
a5 = 0;

x = A(3+n:end,1);
y = A(3+n:end,2);

M2 = [2,1,1;1,2,1;1,1,2];

mata = zeros(npkt,npkt,n);
matb = zeros(npkt,npkt,n);

dborder = A(2,1:end-1);
    
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
    ma = 1/2*M1;
    mb = (d/24)*M2;
    
    for lauf=1:3,
        for lauff=1:3
            %% Beitrag der einzelnen Dreiecke
            mata(ijk(lauff),ijk(lauf),m) = ma(lauff,lauf);
            matb(ijk(lauff),ijk(lauf),m) = mb(lauff,lauf);
        end
    end
    
end

%% Über alle Durchgänge summieren
mata = sum(mata,3);
matb = sum(matb,3);

%% Dirichlet
mata(dborder,:) = [];
mata(:,dborder) = [];

matb(dborder,:) = [];
matb(:,dborder) = [];

[f,lambda] = eigs(-mata,matb,9,'SM');
lambda = -sum(lambda,2);

xnetz = linspace(0,5,101);
ynetz = linspace(0,4,101);
f_inter = zeros(101,101);
c = [1,0,0;-1,1,0;-1,0,1];

for eigen=1:9
    f_eigen = f(:,eigen);
    lambda_eigen = lambda(eigen);
    f_neu = zeros(npkt,1);

    iwert = 1;

    for ind=1:npkt
        if isempty(dborder(dborder==ind)) == 1,
            f_neu(ind) = f_eigen(iwert);
            iwert = iwert + 1;
        end
    end

    f_eigen = f_neu;

    for lauff=1:101
        for lauf=1:101,
            for m=1:n,
                i = A(2+m,1);
                j = A(2+m,2);
                k = A(2+m,3);

                u = ((ynetz(lauff)-y(i))*(x(k)-x(i)) - (y(k)-y(i))*(xnetz(lauf)-x(i)))/((y(j)-y(i))*(x(k)-x(i)) - (y(k)-y(i))*(x(j)-x(i)));
                v = ((ynetz(lauff)-y(i))*(x(j)-x(i)) - (y(j)-y(i))*(xnetz(lauf)-x(i)))/((y(k)-y(i))*(x(j)-x(i)) - (y(j)-y(i))*(x(k)-x(i)));

                if v<=-u+1 && u>=0 && u<=1 && v>=0,
                    cm = c*[f_eigen(i); f_eigen(j); f_eigen(k)];
                    f_inter(lauff,lauf) = cm(1) + cm(2)*u + cm(3)*v;
                end
            end
        end
    end

    subplot(3,3,eigen);
    pcolor(f_inter);
    shading flat;
    axis('equal');
    title(['Eigenschwingung #', num2str(eigen-1)]);
    xlabel('x');
    ylabel('y');
end
