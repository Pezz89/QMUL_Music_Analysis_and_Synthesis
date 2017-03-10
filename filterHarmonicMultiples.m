function [peaks, locs] = filterHarmonicMultiples(peaks, locs)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Remove all peaks that are found to be multiples of the most significant
    % peaks partials.
    %
    % Inputs:
    %  peaks: A matrix of peak frequencies, ordered by most to least
    %  significant
    %  locs: The locations of the peaks in the magnitude spectrum
    % Outputs:
    %  peaks: A filtered matrix of peak frequencies, containing only peaks that
    %  were not calculated to be multiples of others.
    %  locs: The locations of the peaks in the magnitude spectrum
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    if ~isempty(peaks)
        % Create a boolean mask to be used for filtering out peaks
        mask = true(length(peaks), 1);

        % Range in cents
        range = 50;
        % For each peak
        for i=1:length(peaks)
            % If there are no unmasked peaks, skip this iteration
            if ~mask(i)
                continue;
            end

            % Calculate the rounded ratio between the currently selected peaks
            % and all others
            round_rat = round(peaks(i) ./ peaks(i+1:end));
            % Select only those that are higher in frequency than the current
            % peak
            r = find(round_rat >= 1);
            round_rat = round_rat(r);

            % Calculate the absolute difference in cents between each peak and
            % the closest multiple of the current peak
            cent_diff = abs(1200 * log2 (peaks(i+r) ./ peaks(i).*round_rat));

            % Remove cadidates that are within 50 cents of a multiple of the
            % current peak
            removeCandidates = find(cent_diff < range);
            mask(i+removeCandidates) = false;

            i = i + 1;
        end
        % Apply masks to the peaks and their locations
        peaks = peaks(mask);
        locs = locs(mask);
    end
