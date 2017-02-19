
function timeStretch(fileName, ratio)
    % function timeStretch(fileName, ratio)
    %     (based on DAFx Book, ch08/VX_tstretch_real_pv.m)
    %===== this program performs time stretching
    %===== using the FFT-IFFT approach,
    %===== for real ratio, and using
    %===== w1 and w2 windows (analysis and synthesis)
    %===== WLen is the length of the windows
    %===== hopSize and n2: steps (in samples) for the analysis and synthesis

    if (nargin < 2) || (ratio <= 0)
        error('usage: timeStretch(fileName, ratio)');
    end

    %----- user data -----
    n2           = 512; % analysis step [samples]
    hopSize           = round(n2 / ratio); % synthesis step [samples]
    WLen         = 2048; % Window length

    % Read input to be stretched and get relevant meta-data
    [DAFx_in,FS] = audioread(fileName);
    inputInfo = (audioinfo(fileName));
    channels = inputInfo.NumChannels;

    % Sum to mono
    if channels > 1
        DAFx_in = sum(DAFx_in,2);
    end

    % Calculate the length of the input in samples
    L            = length(DAFx_in);
    % Zero-pad input to allow for windowed analysis accross entire file and
    % normalise.
    DAFx_in      = [zeros(WLen, 1); DAFx_in; ...
       zeros(WLen-mod(L,hopSize),1)] / max(abs(DAFx_in));

    [analysis, winCount] = calculateAnalysis(DAFx_in, WLen, hopSize);

    delta = 0.2;
    analysis = normaliseAnalysis(analysis, delta);


    [stable, stable_ratio] = getStable(DAFx_in, analysis, WLen, delta, hopSize);

    timeStretchStable(DAFx_in, FS, stable, ratio / stable_ratio);

function timeStretchStable(in, FS,  stable, ratio)
    %----- time stretching initializations -----
    n2           = 256; % analysis step [samples]
    n1           = round(n2 / ratio); % synthesis step [samples]
    WLen         = 2048; % Window length
    w1           = hanning(WLen); % Hanning window of length WLen
    w2           = w1;
    % TODO; add semi-colon
    tstretch_ratio = n2/n1
    out = zeros(WLen+ceil(length(in)*tstretch_ratio),1);
    omega    = 2*pi*n1*[0:WLen-1]'/WLen;
    phi0     = zeros(WLen,1);
    psi      = zeros(WLen,1);

    devcent = 2*pi*n1/WLen;

    tic
    %UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
    pin  = 0;
    pout = 0;
    pend = length(in)-WLen;

    while pin<pend
        grain = in(pin+1:pin+WLen).* w1;

        if(any(pin+WLen/2 > stable(:, 1) & pin+WLen/2 < stable(:, 2)))
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
            pin  = pin + n1;
            pout = pout + n2;
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
            pin  = pin + n1;
            pout = pout + n1;
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
    system(['open ' outName]);

function [stable, ratio] = getStable(in, analysis, delta, WLen, hopSize)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function to convert transience start and end times to stable part
    % segments.
    % Returns:
    % stable: a 2xN array of segment start and end times, where N is the number
    % of stable parts in the audio.
    % ratio: The ratio between the total size of stable and unstable parts in
    % the audio.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Enables the saving of variables to mat files for plotting in Python
    pythonPlot = true;
    delta = 0

    winCount = floor((length(in)-WLen)/hopSize);
    a = double(analysis > delta)

    % TODO: Convert to Python
    if(false)
        figure
        plot(DAFx_in)
        hold on;
        plot(((1:winCount)*hopSize)+WLen/2,a)
        hold on;
        plot(((1:winCount)*hopSize)+WLen/2,analysis)
    end

    % Code adapted from https://uk.mathworks.com/matlabcentral/newsreader/view_thread/151318
    krn=[1 -1];
    changes=conv(krn, a)
    % Calculate start and end window indexes of transient segments
    t_s = find(changes==1)
    t_e = find(changes==-1);

    % Convert window index to samples
    % TODO: Check sample accuracy of this...
    transience_s = t_s * hopSize;
    transience_e = t_e * hopSize;

    if(pythonPlot)
        % Export variables to mat file for plotting in Python
        s1.analysis = analysis';
        s1.transience_e = transience_e';
        s1.transience_s = transience_s';
        s1.a = a';
        s1.win_count = winCount;
        s1.n1 = hopSize;
        s1.WLen = WLen;
        save('./vars.mat','-struct', 's1')
    end

    % Convert column vectors to rows
    transience_s = transience_s(:);
    transience_e = transience_e(:);

    % Get the length of the input audio
    L = length(in);
    if(transience_s(1) ~= 0)
        transience_e = [0; transience_e];
        transience_s = [transience_s; L];
    end
    stable = horzcat(transience_e, transience_s);
    stable(:, 2) - stable(:, 1);
    ratio = sum(stable(:, 2) - stable(:, 1)) / L;

function [analysis, winCount] = calculateAnalysis(in, WLen, hopSize)
    %----- transience analysis initialization -----
    test = 0.4;
    % Allocate memory to store the current grain to be analysed
    grain = zeros(WLen,1);
    % Allocate memory to store the previous window's magnitude during analysis
    w1           = hanning(WLen); % Analysis Hanning window of length WLen
    mag1 = zeros(WLen/2,1);
    % Calculate the total number of windows to be analysed
    winCount = floor((length(in)-WLen)/hopSize);
    % Allocate memory to store outut analysis values
    analysis = zeros(winCount, 1);

    pin  = 0;
    pout = 1;
    pend = length(in)-WLen;
    %----- transience analysis -----
    while pout<winCount
        grain = in(pin+1:pin+WLen).* w1;
        f = fft(grain);
        mag = abs(f(1:WLen/2));

        analysis(pout) = sqrt(sum((mag-mag1).^2))/(WLen/2);
        %mag_diff = mag-mag1
        %analysis(pout) = sum(mag_diff-abs(mag_diff)/2);
        mag1 = mag;

        pin  = pin + hopSize;
        pout = pout + 1;
    end

function analysis = normaliseAnalysis(analysis, delta)
    % Normalize analysis
    analysis = analysis - mean(analysis);
    analysis = analysis / max(abs(std(analysis)));

    %TODO: Check that this aligns with the original analysis and with the
    %audio.
    analysis = filter(ones(1, 40)/40, 1, analysis);
    analysis = analysis(25:end);

    %TODO: Check that this aligns with the original analysis and with the
    %audio.
    thresh = medfilt1(analysis, 1000);

    analysis = analysis - (delta+thresh);
