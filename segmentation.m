function [V,pitches2]=segmentation(voiced, M, pitches)
% M-file 9.19
% segmentation.m
%
% voiced: original voiced/unvoiced detection
% M: minimum number of blocks in a row with same voiced flag
% pitches: original pitches
% V: changed voiced flag
% pitches2: changed pitches

blocks=length(voiced); % get number of blocks
pitches2=pitches;
V=voiced;
Nv=length(V);

%%%%%%%%%%% step1: eliminate too short voiced segments:
V(Nv+1)=~V(Nv); % one more change at end to get length of last segment
dv=[0, diff(V)]; % derivative
idx=find(dv~=0); % changes in voiced
di=[idx(1)-1,diff(idx)];  % sequence lengths
v0=V(1); % status of 1st sequence
k0=1;   
ii=1; % counter for sequences, idx(ii)-1 is end of sequence
if v0==0
  k0=idx(1); %start of voiced
  ii=ii+1; % first change voiced to unvoiced
end  
while ii<=length(idx);
  L=di(ii);
  k1=idx(ii)-1; % end of voiced sequence
  if L<M
     V(k0:k1)=zeros(1,k1-k0+1);
  end
  if ii<length(idx)
    k0=idx(ii+1); %start of next voiced sequence 
  end
  ii=ii+2;
end

%%%%%%%%%%% step2: eliminate too short unvoiced segments:
V(Nv+1)=~V(Nv); % one more change at end
dv=[0, diff(V)];
idx=find(dv~=0); % changes in voiced
di=[idx(1)-1,diff(idx)];  % sequence lengths
if length(idx)>1 % changes in V
  v0=V(1); % status of 1st sequence
  k0=1;   
  ii=1; % counter for sequences, idx(ii)-1 is end of sequence
  if v0==0
    k0=idx(2); %start of unvoiced
    ii=ii+2; % first change unvoiced to voiced
  end  
  while ii<=length(idx);
    L=di(ii);
    k1=idx(ii)-1; % end of unvoiced sequence
    if L<M
       if k1<blocks % NOT last unvoiced sequence
         V(k0:k1)=ones(1,k1-k0+1);
         % linear pitch interpolation:
         p0=pitches(k0-1);
         p1=pitches(k1+1);
         N=k1-k0+1;
         pitches2(k0:k1)=(1:N)*(p1-p0)/(N+1)+p0;
      end   
    end
    if ii<length(idx)
      k0=idx(ii+1); %start of next unvoiced sequence 
    end
    ii=ii+2;
  end
end  

V=V(1:Nv); % cut last element
%pitches2=pitches2.*V; % set pitches to zero for unvoiced 






