function fourier = fft_ratschek(x)

% M. Ratschek   6-11-2007

% Um auch bei den Matlab-Programmen fft.m und ifft.m dieselbe 
% Reihenfolge der Fourierkoeffizienten zu erhalten wie bei
% den im Skriptum praesentierten C-Programm (four1.c) bzw.
% F90-Programm (fft.f90), hat Herr Ratschek dieses Programm
% geschrieben, bei welchem NACH der Anwendung von fft.m
% "automatisch" eine Umordnung der Matlab-Ergebnisse erfolgt.

% Wichtig ist, dass Sie die Programme fft_ratschek.m und ifft_ratschek.m
% paarweise verwenden, also nicht fft_ratschek.m mit ifft.m
% kombinieren und umgekehrt.

%  n=length(x);
  d=fft(x);

  fourier = cat(find(size(d)>1),d(1),d(end:-1:2));
