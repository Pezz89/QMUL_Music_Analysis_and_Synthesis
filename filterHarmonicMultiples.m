function [peaks, locs] = filterHarmonicMultiples(peaks, locs)
    if ~isempty(peaks)
        mask = true(length(peaks), 1);

        % Range in cents
        range = 50;
        for i=1:length(peaks)
            if ~mask(i)
                continue;
            end

            round_rat = round(peaks(i) ./ peaks(i+1:end));
            r = find(round_rat >= 1);
            round_rat = round_rat(r);

            cent_diff = abs(1200 * log2 (peaks(i+r) ./ peaks(i).*round_rat));

            removeCandidates = find(cent_diff < range);
            mask(i+removeCandidates) = false;

            i = i + 1;
        end
        peaks = peaks(mask);
        locs = locs(mask);
    end
