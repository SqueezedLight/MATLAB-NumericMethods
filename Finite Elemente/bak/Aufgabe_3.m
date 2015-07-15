clear all;
clc;

A = textread('./Rechteck.266');

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

[f,lambda] = eigs(-mata,matb,9,'SM');
lambda = -sum(lambda,2);

xnetz = linspace(0,5,101);
ynetz = linspace(0,4,101);
f_inter = zeros(101,101);
c = [1,0,0;-1,1,0;-1,0,1];
f = f(:,2);
lambda = lambda(2);

for lauff=1:101
    for lauf=1:101,
        for m=1:n,
            i = A(2+m,1);
            j = A(2+m,2);
            k = A(2+m,3);
            
            %v = (xnetz(lauf)*(y(i) - y(j)) + ynetz(lauff)*(x(j) - x(i)) - y(i)*x(j) + x(i)*y(j))/(y(k)*(x(j)-x(i)) + y(i)*(x(k) - x(j)) + y(j)*(x(i) - x(k)));
            %v = (ynetz(lauff)*(x(j) - x(i)) - y(i)*(x(j) - xnetz(lauf)) - y(j)*(xnetz(lauf) - x(i)))/(y(i)*(x(k) - x(j) - 1) + y(j)*(x(i) - x(k)) + y(k)*(x(j) - x(i) +1));
            v = (ynetz(lauff)*(x(j) - x(i)) - y(i)*(x(j) - x(i)) - (y(j) - y(i))*(xnetz(lauf) - x(i)))/((x(j) - x(i))*(y(k) - y(i)) - (y(j) - y(i))*(x(k) - x(i)));
            u = (xnetz(lauf) - x(i) - (x(k) - x(i))*v)/(x(j) - x(i));
            %u = (y(k)*(x(j)*xnetz(lauf) - x(i)*xnetz(lauf) - x(j)*x(i) + x(i)^2) + y(i)*(x(j)*x(k) - x(j)*xnetz(lauf) - x(k)*x(i) + xnetz(lauf)*x(i)) + ynetz(lauff)*(x(i)*x(k) + x(j)*x(i) - x(j)*x(k) - x(i)^2))/(y(i)*(x(k)*x(j) - x(j)^2 - x(k)*x(i) + x(j)*x(i)) + y(j)*(x(i)*x(j) - x(k)*x(j) + x(k)*x(i) - x(i)^2) + y(k)*(x(j) - x(i))^2);
                
            if v<=-u+1 && u>=0 && u<=1 && v>=0,
                cm = c*[f(i); f(j); f(k)];
                f_inter(lauff,lauf) = cm(1) + cm(2)*u + cm(3)*v;
            end
        end
    end
end
f_inter = fliplr(f_inter);    
pcolor(f_inter);
axis('equal');
