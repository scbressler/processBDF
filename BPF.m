function [y,Hd] = BPF(x, fs, norder, cf1, cf2)

% y = randn(1,10*fs)*0.1; % make 10 seconds noise
% n = 10000;
% n = 1000;
% n = 100;
n = norder;

Ny = fs/2;
Wn = [cf1 cf2]/(Ny);
Hd = fir1(n,Wn);

y = flipud(fftfilt(Hd,flipud(fftfilt(Hd,x))));