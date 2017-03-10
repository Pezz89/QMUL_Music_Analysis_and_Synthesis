function [peaks, locs] = filterHarmonicMultiples(peaks, locs)
    if ~isempty(peaks)
        mask = true(length(peaks));
        % Range in cents
        range = 50
        for i=1:length(peaks)
            if ~mask(i)
                continue;
            end

            round_rat = round(peaks(i) / peaks(i+1:end))

            cent_diff = 1200 * log2 (peaks(i+1:end) / peaks(i)*round_rat)
            mask(abs(cent_diff) < range
            keyboard
            
            i = i + 1;
        end
    end
