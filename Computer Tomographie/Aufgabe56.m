%% Definition der Parameter
clear all;
clc;

xi_max = 1.6;

x_min = -1;
y_min = -1;
x_max = 1;
y_max = 1;

n_x = 100;
n_y = 100;

d_x = (x_max - x_min)/n_x;
d_y = (y_max - y_min)/n_y;

x_null = zeros(n_x+1,1);
y_null = zeros(n_y+1,1);
xi_null = zeros(n_y+1);

reply = input('Wellenzahlverschiebung ein? y/n [y]: ', 's');
if isempty(reply),reply = 'y'; end

for rohr=0:3

if rohr == 0, dat = load('./Daten-Rohre/rohre.0'); end
if rohr == 1, dat = load('./Daten-Rohre/rohre.1'); end
if rohr == 2, dat = load('./Daten-Rohre/rohre.2'); end
if rohr == 3, dat = load('./Daten-Rohre/rohre.3'); end

nphi = length(dat(:,1));
n = length(dat(1,:));
delta_xi = 2*xi_max/n;
p_int = zeros(n_y+1,nphi);

%% Ausserste Schleife ueber die Winkel PHI
for i=1:nphi
    
phi = i*(pi/nphi);
k_xi = zeros(n,1);
xi = zeros(n,1);
    
%% Schleife ueber die verschiedenen xi Positionen
for j=1:n       
     xi(j) = -xi_max + (j-1)*delta_xi; 

        if reply == 'y'       
            if (1<=j  && j<(n/2)+1), k_xi(j) = 2*pi*(j-1)/(n*delta_xi); end
            if ((n/2)+1<=j && j<=n), k_xi(j) = 2*pi*(j-1-n)/(n*delta_xi); end      
        end
     
        if reply ~= 'y'       
            k_xi(j) = 2*pi*(j-1)/(n*delta_xi); 
        end       
end

p = dat(i,:);
p = p';



%% Fourier Koeffizienten berechnen
fkoeff = fft_ratschek(p)*(n/32);
fkoeff = fkoeff.*abs(k_xi);
ifourier = ifft_ratschek(fkoeff);   

%% Ergebnismatrix
for ii=1:n_x + 1
    x_null(ii) = x_min + (ii - 1)*d_x;
    y_null(ii) = y_min + (ii - 1)*d_y;
end

for rr=1:n_x+1
    for ss=1:n_x+1
        xi_null(rr,ss) = x_null(ss) * cos(phi) + y_null(rr) * sin(phi);
        ind = find(xi_null(rr,ss) > xi);
        ind = ind(end);

        p_int(rr,ss,i) = 10 * ((xi_null(rr,ss) - xi(ind))*ifourier(ind + 1) + (xi(ind + 1) - xi_null(rr,ss))*ifourier(ind));
    end
end

end

fimage = 1/(2*nphi)*sum(p_int,3);

if rohr == 0, h = figure; end
subplot(2,2,rohr+1);
contour(fimage,25);
title(['Contour-Plot von rohr.', num2str(rohr), ' mit n= ', num2str(n), ' nphi=' num2str(nphi)]);
xlabel('x - Index');
ylabel('y - Index');
axis('square');

end
