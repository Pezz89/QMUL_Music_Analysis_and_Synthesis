function [f0 B] = fundamentals(audioFile, transcriptionMatrix)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This program estimates the fundamental frequency of notes in a
    % harpsichord recording, given an audio file path as input.
    % Based on transcription approach in  (Dixon, Mauch and Tidhar, 2011)
    % "Estimation of harpsichord inharmonicity and temperament from musical
    % recordings"
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
            % Get the current peak's index in the spectrum
            binInd = locs(i);
            % Calculate variables for quadratic interpolation
            ap1 = X(binInd+1);
            am1 = X(binInd-1);
            a = X(binInd);
            % Interpolate peaks to refine peak frequency
            delta = 0.5*(am1 - ap1)/(2*(am1 - 2*a + ap1));
            peaks(i) = (binInd+delta) * binFreq;
            i = i + 1;
        end
        % Apply "Conservative transcription" harmonic filtering to remove
        % multiples of peaks assumed to be f0s
        [peaks, locs] = filterHarmonicMultiples(peaks, locs);

        % Store peaks with previously calculated peaks
        allPeaks = [allPeaks, peaks];
        % Increment read pointer
        pin = pin + hop;
    end


    % Define nominal pitch
    fst = 415;
    % Calculate distance in cents between the nominal pitch and all found peaks
    cent_diff = abs(1200 * log2 (allPeaks ./ fst));
    % Calculate the average frequency of peaks that are within a semi tone of
    % the nominal pitch
    fst = median(allPeaks(cent_diff < 50));
    % Generate a range of pitches for midi notes based on the calculated
    % nominal pitch. This equation was adapted from: https://www.midikits.net/midi_analyser/midi_note_frequency.htm<Paste>
    % Analysed notes range from MIDI note 36 (C2) to 80 (G#5) as per the Dixon
    % et. al paper
    pitchRange = -33:11;
    f = 2.^(pitchRange/12).*fst;

    % Find all peaks that are within a semi tone of each midi note, and
    % calculate their average to generate the estimation of frequency
    i = 1;
    % Allocate empty matrix for each note
    freq = zeros(length(f), 1);
    while i < length(f)
        cent_diff = abs(1200 * log2 (allPeaks / f(i)));
        % If no peaks can be found within the range, the note is zerod
        if ~any(cent_diff)
            freq(i) = 0;
        end
        % Calculate the median, as above...
        freq(i) = median(allPeaks(cent_diff < 50));
        i = i + 1;
    end
    % Fill all notes that are not analysed with zeros to create the final
    % matrix
    f0 = [zeros(35, 1); freq; zeros(47, 1)];

    % Due to time constraints, it was not possible to implement the
    % inharmonicity values
    B = [];
end
