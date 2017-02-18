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
    win_count = floor((length(DAFx_in)-WLen)/n1)
    analysis = zeros(win_count, 1);

    pin  = 0;
    pout = 1;
    pend = length(DAFx_in)-WLen;
    %----- transience analysis -----
    while pout<win_count
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
    analysis = analysis / max(analysis);

    % TODO: Absolute values seems odd and possibly wrong... check this...
    thresh = zeros(length(analysis)-2, 1);
    for i = 3:length(analysis)
        thresh(i-2) = mean( ...
            [abs(analysis(i)), ...
            abs(analysis(i-1)),...
            abs(analysis(i-2))]...
        );
    end
    delta = 0.1;
    lambda = 0.0;

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

    % Export variables to mat file for plotting in Python
    win_count
    s1.analysis = analysis';
    s1.transience_e = transience_e' ;
    s1.transience_s = transience_s';
    s1.a = a';
    s1.win_count = win_count;
    s1.n1 = n1;
    s1.WLen = WLen;
    save('./vars.mat','-struct', 's1')

    [stable, stable_ratio] = getStable(DAFx_in, transience_s, transience_e);

    timeStretchStable(DAFx_in, FS, stable, ratio / stable_ratio);

function timeStretchStable(in, FS,  stable, ratio)
    %----- time stretching initializations -----
    n2           = 512; % analysis step [samples]
    n1           = round(n2 / ratio); % synthesis step [samples]
    WLen         = 2048; % Window length
    w1           = hanning(WLen); % Hanning window of length WLen
    w2           = w1;
    tstretch_ratio = n2/n1
    out = zeros(WLen+ceil(length(in)*tstretch_ratio),1);
    omega    = 2*pi*n1*[0:WLen-1]'/WLen;
    phi0     = zeros(WLen,1);
    psi      = zeros(WLen,1);

    devcent = 2*pi*n1/WLen;

    tic
    %UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
    pin  = -n1;
    pout = -n1;
    pend = length(in)-WLen;

    while pin<pend
        if(any(pin >= stable(:, 1) & pin+WLen < stable(:, 2)))
            pin  = pin + n1;
            pout = pout + n2;
        else
            pin  = pin + n1;
            pout = pout + n1;
        end
        if(pin>=pend)
            break;
        end

        grain = in(pin+1:pin+WLen).* w1;

        if(any(pin >= stable(:, 1) & pin+WLen < stable(:, 2)))
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
            out(pout+1:pout+WLen) = ...
               out(pout+1:pout+WLen) + grain;

        else
            f     = fft(fftshift(grain));
            r     = abs(f);
            phi   = angle(f);
            delta_phi= omega + princarg(phi-phi0-omega);
            phi0  = phi;
            psi   = princarg(psi+delta_phi);
            ft    = (r.* exp(i*psi));
            grain = fftshift(real(ifft(ft))).*w2;
            out(pout+1:pout+WLen) = ...
               out(pout+1:pout+WLen) + grain/tstretch_ratio;
        end
    end
    %UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
    toc

    %----- listening and saving the output -----
    %in  = in(WLen+1:WLen+L);
    out = out(WLen+1:length(out))/max(abs(out));
    % soundsc(out, FS);
    outName = ['./out' sprintf('%3.1f', ratio) '.wav'];
    wavwrite(out, FS, outName);
    system(['play --silent ' outName]);

function [stable, ratio] = getStable(in, transience_s, transience_e)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function to convert transience start and end times to stable part
    % segments.
    % Returns:
    % stable: a 2xN array of segment start and end times, where N is the number
    % of stable parts in the audio.
    % ratio: The ratio between the total size of stable and unstable parts in
    % the audio.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Convert column vectors to rows
    transience_s = transience_s(:);
    transience_e = transience_e(:);

    % Get the length of the input audio
    L = length(in);
    if(transience_s(1) ~= 0)
        transience_e = [0; transience_e];
        transience_s = [transience_s; L];
    end
    stable = horzcat(transience_e, transience_s)
    stable(:, 2) - stable(:, 1)
    ratio = sum(stable(:, 2) - stable(:, 1)) / L
