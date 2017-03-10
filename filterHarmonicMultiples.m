function [peaks, locs] = filterHarmonicMultiples(peaks, locs)
    if ~isempty(peaks)
        mask = true(length(peaks), 1);
        [peaks, inds] = sort(peaks);
        peaks = flip(peaks);
        locs = locs(inds);
        locs = flip(locs);

        % Range in cents
        range = 50;
        for i=1:length(peaks)
            if ~mask(i)
                continue;
            end

            round_rat = round(peaks(i) ./ peaks(i+1:end));

            cent_diff = abs(1200 * log2 (peaks(i+1:end) ./ peaks(i).*round_rat));

            removeCandidates = find(cent_diff < range);
            mask(i+removeCandidates) = false;

            i = i + 1;
        end
        peaks = peaks(mask);
        locs = locs(mask);
    end
