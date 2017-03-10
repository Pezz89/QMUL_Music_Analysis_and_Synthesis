function [peaks, locs] = filterHarmonicMultiples(peaks, locs)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Remove all peaks that are found to be multiples of the most significant
    % peaks partials.
    %
    % Inputs:
    %  peaks: A matrix of peak frequencies
    %  locs: The locations of the peaks in the magnitude spectrum
    % Outputs:
    %  peaks: A filtered matrix of peak frequencies, containing only peaks that
    %  were not calculated to be multiples of others.
    %  locs: The locations of the peaks in the magnitude spectrum
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
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
