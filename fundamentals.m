function [f0 B] = fundamentals(audioFile, transcriptionMatrix)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Inputs:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % analysis step [samples]
    n1 = 256;
    % synthesis step [samples]
    n2 = n1;
    % Window length
    WLen = 2048;
    % Hanning window of length WLen
    w1 = hanning(WLen);

    [in, fs] = wavread(audioFile);
    %----- initialize windows, arrays, etc -----
    window = hanning(WLen, 'periodic'); % input window
    nChannel = WLen/2;
    L = length(in);
    % Convert to mono
    in = (in(:, 1)*.5)+(in(:, 2)*.5);
    in = [zeros(WLen, 1); in; ...
    zeros(WLen-mod(L,n1),1)] / max(abs(in));

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
