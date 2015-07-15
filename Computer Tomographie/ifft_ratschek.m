function ifourier = ifft_ratschek(d)

% M. Ratschek   6-11-2007

% Um auch bei den Matlab-Programmen fft.m und ifft.m dieselbe 
% Reihenfolge der Fourierkoeffizienten zu erhalten wie bei
% den im Skriptum praesentierten C-Programm (four1.c) bzw.
% F90-Programm (fft.f90), hat Herr Ratschek dieses Programm
% geschrieben, bei welchem VOR der Anwendung von ifft.m
% "automatisch" eine Umordnung der Matlab-Ergebnisse erfolgt.

% Wichtig ist, dass Sie die Programme ifft_ratschek.m und fft_ratschek.m
% paarweise verwenden, also nicht ifft_ratschek.m mit fft.m
% kombinieren und umgekehrt.

%  n=length(d);
  x=cat(find(size(d)>1),d(1),d(end:-1:2));

  ifourier=ifft(x);
