function [H, phi]  = fftdb(x, Nfft)
% [H, phi]  = fftdb(x, Nfft)
%
% M-file 9.16
% fftdb.m
%
% compute FFT in db
% x: input signal
% Nfft: FFT length
% H: FFT magnitude values in dB
% phi: phase values
% 
% (c) 2002 Florian Keiler

if nargin<2
  Nfft = length(x);
end

F      = fft(x,Nfft);
F      = F(1:Nfft/2+1);       % f=0,..,Fs/2
phi    = angle(F);            % phase in [-pi,pi]
F      = abs(F)/Nfft*2;       % normalize to FFT length
%----- return -100 db for F==0 to avoid "log of zero" warnings -----
H      = -100*ones(size(F));
idx    = find(F~=0);
H(idx) = 20*log10(F(idx));    % non-zero values in dB
