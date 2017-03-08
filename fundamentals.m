function [f0 B] = fundamentals(audioFile, transcriptionMatrix)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Inputs:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    stretchRatio = stretchRatio / stableRatio;
    % analysis step [samples]
    n2 = 256;
    % synthesis step [samples]
    n1 = round(n2 / stretchRatio);
    % Window length
    WLen = 2048;
    % Hanning window of length WLen
    w1 = hanning(WLen);
    % Calculate ratio between analysis and synthesis hop size as the stretch
    % ratio
    tstretch_ratio = n2/n1;

    % Allocate memory for output samples
    out = zeros(WLen+ceil(length(in)*(1-stableRatio)) + WLen*2+ceil(length(in)*stableRatio*stretchRatio),1);
    length(out)
    % Initialize memory for phase vocoder variables
    omega    = 2*pi*n1*[0:WLen-1]'/WLen;
    phi0     = zeros(WLen,1);
    psi      = zeros(WLen,1);

    % Initialize read and write pointers for audio input and output
    pin  = 0;
    pout = 0;
    % Calculate the length of input audio
    pend = length(in)-WLen;

    while pin<pend
        % Read grain from input and apply hanning window
        grain = in(pin+1:pin+WLen).* w1;

        %===========================================
        f     = fft(fftshift(grain));
        r     = abs(f);
        phi   = angle(f);
        delta_phi= omega + princarg(phi-phi0-omega);
        phi0  = phi;
        fi = (1/(2*pi))*(delta_phi/n1)*fs;
        % ===========================================
        % Increament read and write pointers by hope sizes
        pin  = pin + n1;
        pout = pout + n2;
    end

    % Normalise output
    out = out(WLen+1:length(out))/max(abs(out));
    % Write audio out and open in the deafult system application
    outName = ['./out' sprintf('%3.1f', stretchRatio) '.wav'];
    wavwrite(out, FS, outName);
    %system(['open ' outName]);
