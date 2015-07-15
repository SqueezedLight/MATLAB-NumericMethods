%% Definition der Parameter
clear all;
clc;

nphi = [12,24,48,96];
n = 32;
xi_max = 1.6;
delta_xi = 2*xi_max/n;

p = zeros(n,1);

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

%% Ausserste Schleife ueber die Winkel PHI

for pp=1:numel(nphi)

p_int = zeros(n_y+1,nphi(pp));  
    
for i=1:nphi(pp)
    
    phi = i*(pi/nphi(pp));
    k_xi = zeros(n,1);
    xi = zeros(n,1);
    
%% Schleife ueber die verschiedenen xi Positionen
   for j=1:n
       
        xi(j) = -xi_max + (j-1)*delta_xi; 
        p(j) = proj_test1(xi(j),phi);

        if reply == 'y'       
            if (1<=j  && j<(n/2)+1), k_xi(j) = 2*pi*(j-1)/(n*delta_xi); end
            if ((n/2)+1<=j && j<=n), k_xi(j) = 2*pi*(j-1-n)/(n*delta_xi); end      
        end
     
        if reply ~= 'y'       
            k_xi(j) = 2*pi*(j-1)/(n*delta_xi); 
        end
   end

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

fimage = 1/(2*nphi(pp))*sum(p_int,3);
if pp ==1, g = figure; end
subplot(2,2,pp);
contour(fimage,25);
title(['Contour - Plot: n=32   nphi= ', num2str(nphi(pp))]);
xlabel('x - Index');
ylabel('y - Index');
axis('square');

end