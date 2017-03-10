
bach_pr1;

%Load audiofile
[Y, FS] = audioread('bach_pr1_a.wav');
% Sum to mono
Y = sum(Y,size(Y,2))*0.5;

%Set window parameters
Wlen = 4096;
hop = 512;
w = blackmanharris(Wlen);
binFreq = (FS/2) / Wlen;

%Sample at 10 second mark
pin = 0;
pend = length(Y)-Wlen-1;

midiNoteFreqs = zeros(127, 1);
midiNoteCounts = zeros(127, 1);
while pin < pend
    grain = Y(pin+1:pin+Wlen) .* w;

    f = fft(grain);
    X = abs(f);
    X = X(1:length(X)/2);
    [peaks, locs] = findPeaks(X, binFreq);
    peaks = peaks';
    [peaks, locs] = filterHarmonicMultiples(peaks, locs);
    peaks = peaks(1:min(length(peaks), 40));
    locs = locs(1:min(length(peaks), 40));

    midPoint = (pin+Wlen/2)/FS;

    notesInRange = midPoint > notes(:,1) & midPoint < notes(:,2);
    currentNotes = notes(notesInRange,3);

    % Get midi note frequency
    f = (440/32)*2.^((currentNotes-9)/12);

    if ~isempty(f) & length(peaks) >= length(f)
        [M, I] = min(abs(f - peaks)');
        intFreq = zeros(length(f), 1);
        i = 1;
        while i < length(f)
            binInd = locs(I(i));
            ap1 = X(binInd+1);
            am1 = X(binInd-1);
            a = X(binInd);
            delta = 0.5*(am1 - ap1)/(2*(am1 - 2*a + ap1));
            intFreq(i) = (binInd+delta) * binFreq;
            i = i + 1;
        end
        if ~isempty(intFreq)
            midiNoteFreqs(currentNotes) = midiNoteFreqs(currentNotes) + intFreq;
            midiNoteCounts(currentNotes) = midiNoteCounts(currentNotes) + 1;
        end
        %scatter(ones(length(intFreq), 1) * (pin * FS), intFreq, 'filled');
        %hold on;
    end

    %axis tight;
    pin = pin + hop;
end

midiNoteFreqs(midiNoteFreqs ~= 0)  = midiNoteFreqs(midiNoteFreqs ~= 0) ./ midiNoteCounts(midiNoteFreqs ~= 0)

