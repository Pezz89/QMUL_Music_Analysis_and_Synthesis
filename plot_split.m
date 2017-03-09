function plot_slit(idx, t, x)
% plot_slit(idx, t, x)
% 
% M-file 9.20
% plot_split.m
%
% plot segments of a vector
% idx: vector with positions of vector x to be plotted
% t: time vector
% x: signal vector
% x is segmented into parts
% 
% (c) 2002 Florian Keiler

di=diff(idx);
L=length(di);

n0=1;
pos_di=find(di>1);
ii=1; % counter for pos_di

hold off
while ii<=length(pos_di) %n0<=length(x)
  n1=pos_di(ii);
  plot(t(idx(n0:n1)),x(idx(n0:n1)))
  hold on
  n0=n1+1;
  ii=ii+1;
end  

n1=length(idx);
plot(t(idx(n0:n1)),x(idx(n0:n1)))
hold off