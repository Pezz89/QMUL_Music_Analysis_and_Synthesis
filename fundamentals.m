function [f0 B] = fundamentals(audioFile, transcriptionMatrix)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This program estimates the fundamental frequency of notes in a
    % harpsichord recording, given an audio file path as input.
    %
    % Inputs:
    %   audioFile: Input harpsichord audio file path
    %   transcriptionMatrix: Not used
    % Outputs:
    %  f0: A matrix of 127 values, representing the fundamental frequency of
    %  each midi note. Unused values are set to 0.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Load audiofile
    [Y, FS] = audioread(audioFile);
    % Sum to mono
    Y = sum(Y,size(Y,2))*0.5;

    % Set window parameters
    Wlen = 4096;
    % Set hop size
    hop = 512;
    % Generate window function
    w = blackmanharris(Wlen);
    % Calculate the frequency incremenet between consecutive FFT bins
    binFreq = (FS/2) / Wlen;

    % Create start and end pointers for input audio
    pin = 0;
    pend = length(Y)-Wlen-1;

    % Create empty arrays to store f0 peaks
    allPeaks = [];
    while pin < pend
        % Get a single grain and apply windowing function
        grain = Y(pin+1:pin+Wlen) .* w;

        % Generate FFT of the grain for spectral analysis
        X = fft(grain);
        % Calculate the magnitude
        X = abs(X);
        % Take only the non-mirrored FFT bins
        X = X(1:length(X)/2);

        % Get ordered matrix of peaks that represent harmonic content
        [peaks, locs] = findPeaks(X, binFreq);
        peaks = peaks';

        % Get the log magnitude
        X = log(X);

        i = 1;
        % For each peak found...
        while i < length(peaks)

            binInd = locs(i);
            ap1 = X(binInd+1);
            am1 = X(binInd-1);
            a = X(binInd);
            delta = 0.5*(am1 - ap1)/(2*(am1 - 2*a + ap1));
            peaks(i) = (binInd+delta) * binFreq;
            i = i + 1;
        end
        % "Conservative transcription" harmonic filtering.
        [peaks, locs] = filterHarmonicMultiples(peaks, locs);

        allPeaks = [allPeaks, peaks];
        pin = pin + hop;
    end


    % Define nominal pitch
    fst = 415;
    cent_diff = abs(1200 * log2 (allPeaks ./ fst));
    fst = median(allPeaks(cent_diff < 50));
    pitchRange = -33:11;
    f = 2.^(pitchRange/12).*fst;

    i = 1;
    freq = zeros(length(f), 1);
    while i < length(f)
        cent_diff = abs(1200 * log2 (allPeaks / f(i)));
        if ~any(cent_diff)
            freq(i) = 0;
        end
        freq(i) = median(allPeaks(cent_diff < 50));
        i = i + 1;
    end
    f0 = [zeros(35, 1); freq; zeros(47, 1)];

    B = [];
end
