function [idx, idx0] = find_loc_max(x)
% [idx, idx0] = find_loc_max(x)
%
% M-file 9.17
% find_loc_max.m
%
% find local maxima in vector x
%      idx : positions of local max.
%      idx0: positions of local max. with 2 identical values
%      if only 1 return value: positions of all maxima
% 
% (c) 2002 Florian Keiler

N    = length(x);
dx   = diff(x);            % derivation, we need to find sign changes from + to -
dx1  = dx(2:N-1);
dx2  = dx(1:N-2);
prod = dx1.*dx2;
idx1 = find(prod<0);       % sign change in dx1
idx2 = find(dx1(idx1)<0);  % only change from + to -
idx  = idx1(idx2)+1;       % positions of single maxima
%----- zeros in dx? => maxima with 2 identical values -----
idx3 = find(dx==0);
idx4 = find(x(idx3)>0);    % only maxima
idx0 = idx3(idx4);
%----- positions of double maxima, same values at idx3(idx4)+1 -----
if nargout==1              % output 1 vector with positions of all maxima
   idx = sort([idx,idx0]); % (for double max. only 1st position)
end
