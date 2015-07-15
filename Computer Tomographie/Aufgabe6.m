%% Definition der Parameter
clear all;
clc;

dat = load('./Daten-Rohre/rohre.3');
nphi = length(dat(:,1));
n = length(dat(1,:));
% nphi = 12;
% n = 32;
xi_max = 1.6;
delta_xi = 2*xi_max/n;

p = zeros(1,n);
% for iphi=1:nphi
%     p = dat(iphi,:);
%     p = p';
% end

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
p_int = zeros(n_y+1,nphi);
xi_null = zeros(n_y+1);

reply = input('Wellenzahlverschiebung ein? y/n [y]: ', 's');
if isempty(reply),reply = 'y'; end

%% Ausserste Schleife ueber die Winkel PHI
for i=1:nphi
    
phi = i*(pi/nphi);
k_xi = zeros(n,1);
xi = zeros(n,1);
    
%% Schleife ueber die verschiedenen xi Positionen
for j=1:n       
     xi(j) = -xi_max + (j-1)*delta_xi; 
     %p(j) = proj_test1(xi(j),phi);
        
     if (1<=j  && j<(n/2)+1), k_xi(j) = 2*pi*(j-1)/(n*delta_xi); end
     if ((n/2)+1<=j && j<=n), k_xi(j) = 2*pi*(j-1-n)/(n*delta_xi); end        
end



% for iphi=1:nphi
p = dat(i,:);
p = p';
% end



%% Fourier Koeffizienten berechnen
fkoeff = fft_ratschek(p)*(n/32);
if reply == 'y', fkoeff = fkoeff.*abs(k_xi); end
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
g = figure;
surf(fimage');
%axis([0 n_x 0 n_y -0.3 1]);

h = figure;
contour(fimage',25);
axis('square');