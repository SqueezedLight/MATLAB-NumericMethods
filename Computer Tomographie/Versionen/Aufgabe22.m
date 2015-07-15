%% Definition der Parameter
clear all;
clc;

nphi = 48;
n = [32,64,128,256];
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

%% Ausserste Schleife ueber die Winkel PHI

for pp=1:numel(n)

delta_xi = 2*xi_max/n(pp);    
p_int = zeros(n_y+1,nphi);  
p = zeros(n(pp),1);
    
for i=1:nphi
    
    phi = i*(pi/nphi);
    %phi = phi*180/pi;
    k_xi = zeros(n(pp),1);
    xi = zeros(n(pp),1);
%% Schleife ueber die verschiedenen xi Positionen
   for j=1:n(pp)
       
        xi(j) = -xi_max + (j-1)*delta_xi; 
        p(j) = proj_test1(xi(j),phi);
        
        if (1<=j  && j<(n(pp)/2)+1), k_xi(j) = 2*pi*(j-1)/(n(pp)*delta_xi); end
        if ((n(pp)/2)+1<=j && j<=n(pp)), k_xi(j) = 2*pi*(j-1-n(pp))/(n(pp)*delta_xi); end
        
   end

%% Fourier Koeffizienten berechnen
fkoeff = fft_ratschek(p);  
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
if pp ==1, g = figure; end
subplot(2,2,pp);
% contour(fimage',30);
% axis('square');

% h = figure;
surf(fimage');
axis([0 n_x 0 n_y -0.3 1]);
end