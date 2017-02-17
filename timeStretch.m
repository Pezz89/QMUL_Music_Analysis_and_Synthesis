function timeStretch(fileName, ratio)
% function timeStretch(fileName, ratio)
%     (based on DAFx Book, ch08/VX_tstretch_real_pv.m)
%===== this program performs time stretching 
%===== using the FFT-IFFT approach, 
%===== for real ratio, and using
%===== w1 and w2 windows (analysis and synthesis)
%===== WLen is the length of the windows
%===== n1 and n2: steps (in samples) for the analysis and synthesis

if (nargin < 2) || (ratio <= 0)
	error('usage: timeStretch(fileName, ratio)');
end

%----- user data -----
n2           = 512; % analysis step [samples]
n1           = round(n2 / ratio); % synthesis step [samples]
WLen         = 2048; % Window length
w1           = hanning(WLen); % Hanning window of length WLen
w2           = w1;
[DAFx_in,FS,channels] = wavread(fileName);
if channels > 1
	DAFx_in = sum(DAFx_in,2);
end
L            = length(DAFx_in);
DAFx_in      = [zeros(WLen, 1); DAFx_in; ...
   zeros(WLen-mod(L,n1),1)] / max(abs(DAFx_in));

%----- transience analysis initialization -----
test = 0.4;
devcent = 2*pi*n1/WLen;
vtest = test * devcent;
grain = zeros(WLen,1);
theta1 = zeros(WLen,1);
theta2 = zeros(WLen,1);
mag1 = zeros(WLen/2,1);
mag2 = zeros(WLen/2,1);
win_count = floor((length(DAFx_in)-WLen)/n1);
analysis = zeros(win_count, 1);

pin  = 0;
pout = 1;
pend = length(DAFx_in)-WLen;
%----- transience analysis -----
while pin<pend
    grain = DAFx_in(pin+1:pin+WLen).* w1;
    f = fft(grain);
    mag = abs(f(1:WLen/2));


    analysis(pout) = sqrt(sum((mag-mag1).^2))/(WLen/2);
    %mag_diff = mag-mag1
    %analysis(pout) = sum(mag_diff-abs(mag_diff)/2);
    mag1 = mag;

	pin  = pin + n1;
    pout = pout + 1;

end

% Normalize analysis
analysis = analysis - mean(analysis);
analysis = analysis / max(analysis)

% TODO: Absolute values seems odd and possibly wrong... check this...
thresh = zeros(length(analysis)-2, 1);
for i = 3:length(analysis)
    thresh(i-2) = mean( ...
        [abs(analysis(i)), ...
        abs(analysis(i-1)),...
        abs(analysis(i-2))]...
    );
end
thresh
delta = 0.05;
lambda = 1.00;

a = zeros(length(thresh),1);
a(i) = analysis(1) > thresh(1);
for i = 3:length(thresh)
    if a(i-2) == true
        a(i-1) = analysis(i) > delta - lambda * thresh(i-2);
    else
        a(i-1) = analysis(i) > delta + lambda * thresh(i-2);
    end
end

if(false)
    figure
    %plot(DAFx_in)
    %hold on;
    plot(((1:win_count)*n1)+WLen/2,a)
    hold on;
    plot(((1:win_count)*n1)+WLen/2,analysis)
end

% Code adapted from https://uk.mathworks.com/matlabcentral/newsreader/view_thread/151318
krn=[1 -1];
changes=conv(krn, a);
% Calculate start and end window indexes of transient segments
t_s = find(changes==1);
t_e = find(changes==-1);

% Convert window index to samples
% TODO: Check sample accuracy of this...
transience_s = t_s * n1;
transience_e = t_e * n1;

s1.analysis = analysis';
s1.transience_e = transience_e' ;
s1.transience_s = transience_s';
s1.a = a';
s1.win_count = win_count;
s1.n1 = n1;
s1.WLen = WLen;

save('./vars.mat','-struct', 's1')

if(false)
    figure
    %plot(DAFx_in)
    plot(((1:win_count)*n1)+WLen/2,analysis)
    plot()
end
return

%-------------------------------

%----- time stretching initializations -----
tstretch_ratio = n2/n1
DAFx_out = zeros(WLen+ceil(length(DAFx_in)*tstretch_ratio),1);
omega    = 2*pi*n1*[0:WLen-1]'/WLen;
phi0     = zeros(WLen,1);
psi      = zeros(WLen,1);

devcent = 2*pi*n1/WLen;

tic
%UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
pin  = 0;
pout = 0;
pend = length(DAFx_in)-WLen;
while pin<pend
	grain = DAFx_in(pin+1:pin+WLen).* w1;
%===========================================
	f     = fft(fftshift(grain));
	r     = abs(f);
	phi   = angle(f);
	delta_phi= omega + princarg(phi-phi0-omega);
	phi0  = phi;
	psi   = princarg(psi+delta_phi*tstretch_ratio);
	ft    = (r.* exp(i*psi));
	grain = fftshift(real(ifft(ft))).*w2;
	% plot(grain);drawnow;
% ===========================================
	DAFx_out(pout+1:pout+WLen) = ...
	   DAFx_out(pout+1:pout+WLen) + grain;
	pin  = pin + n1;
	pout = pout + n2;
end
%UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
toc

%----- listening and saving the output -----
%DAFx_in  = DAFx_in(WLen+1:WLen+L);
DAFx_out = DAFx_out(WLen+1:length(DAFx_out))/max(abs(DAFx_out));
% soundsc(DAFx_out, FS);
outName = [fileName(1:end-4) sprintf('%3.1f', ratio) '.wav'];
wavwrite(DAFx_out, FS, outName);
system(['play --silent ' outName]);
