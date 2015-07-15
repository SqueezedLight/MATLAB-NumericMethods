clear all;
clc;

A = textread('./dreiecke_7.txt');

n = A(1,1);
npkt = A(1,2);
inhom = zeros(npkt,1);
mat = zeros(npkt,npkt,n);
mata = zeros(npkt,npkt,n);
matb = zeros(npkt,npkt,n);
cborder = [1,6,2]; 
dborder = [3,4,8];
centry = 1;

h = 3;
a1 = 2;
a2 = 2;
a4 = -0.5;
a5 = 0;

for m=1:n,
    
    i = A(2+m,1);
    j = A(2+m,2);
    k = A(2+m,3);
    ijk = [i,j,k];
    
    x = A(3+n:end,1);
    y = A(3+n:end,2);
    
    d = (x(j)-x(i))*(y(k)-y(i)) - (x(k)-x(i))*(y(j)-y(i));      %Jakobi Determinante
    
    e11 = (a1*(y(k)-y(i))^2 + a2*(x(k)-x(i))^2)/d;
    e12 = -(a1*(y(k)-y(i))*(y(j)-y(i)) + a2*(x(k)-x(i))*(x(j)-x(i)))/d;
    e22 = (a1*(y(j)-y(i))^2 + a2*(x(j)-x(i))^2)/d;
    
    M1 = [e11+2*e12+e22,-e11-e12,-e12-e22;-e11-e12,e11,e12;-e12-e22,e12,e22];
    b = h*d/6.*[1;1;1];
    ma = 1/2*M1;
    
    %% Beitrag der einzelnen Dreiecke
    for lauf=1:3,
        for lauff=1:3
            mata(ijk(lauff),ijk(lauf),m) = ma(lauff,lauf);
            inhom(ijk(lauff),m) = b(lauf);

        end
    end
 
    
    %% Beitrag aus den Causchy'schen Randbedingungen
    if isempty(cborder(cborder==i)) == 0,
        d_pq = sqrt((x(i)-x(j))^2 + (y(i)-y(j))^2);
        mb = a4*d_pq/6.*[2,1;1,2];
        
        for lauf_1=1:2,
            for lauff_1=1:2
                matb(ijk(lauff_1),ijk(lauf_1),m) = mb(lauff_1,lauf_1);
            end
        end
        centry = centry + 1;
    end
    %% Summieren
    mat(:,:,m) = mata(:,:,m) + matb(:,:,m);
    
    
end

mat = sum(mat,3);
inhom = sum(inhom,2);

%% Dirichlet berücksichtigen

wert = 20;

for ind=dborder,
    for np=1:npkt,
        inhom(np) = inhom(np) - mat(np,ind)*wert;
        mat(ind,np) = 0;
        mat(np,ind) = 0;
    end


mat(ind,ind) = 1;
inhom(ind) = wert;

end

%% Temperatur berechnen
T = mat\inhom;

