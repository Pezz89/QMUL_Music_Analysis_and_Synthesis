function [freqs, inds] = findPeaks(X, frequencyPerBin)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find significant harmonic peaks in the magnitude spectrum.
    %
    % Inputs:
    %  X: The magnitude spectrum of a grain
    %  frequencyPerBin: The frequency increment between consecutive FFT bins
    % Outputs:
    %  freqs: A matrix of frequencies for each found peak
    %  inds: FFT indexes of found peaks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find the maximum value in the spectrum
    Xmax = max(X);
    % Create an adaptive threhold using a moving mean and moving standard
    % deviation filter. Also create a static threshold at -25dB below the
    % maximum value in the spectrum.
    ThreshInds = find(X > (movmean(X, 200) + (0.5 * movstd(X, 200))) & X > (10^-2.5)*Xmax);

    % Visualise the threshold for the current grain
    if false
        figure
        plot((movmean(X, 200) + (0.5 * movstd(X, 200))), 'g')
        hold on
        plot(X)
        hold on
    end

    % Check that there are values above the threshold
    if any(ThreshInds)
        % Find all peaks in the spectrum
        [peaks, locs] = findpeaks(X);
        % Get the index of peaks that are above the threshold
        [~, ~, ib] = intersect(ThreshInds, locs);

        % Filter out all peaks that aren't above the threshold
        peaks = peaks(ib);
        locs = locs(ib);
        % Calculate the frequencies for all selected peaks 
        freqs = locs * frequencyPerBin;
        ampfreq = [peaks freqs locs];

        % Sort peaks from lowest to highest amplitude
        sorted = sortrows(ampfreq,1);
        freqs = flip(sorted(:,2));
        inds = flip(sorted(:,3));
        % Filter out peaks at the start or end as they cannot be interpolated
        freqs = freqs(inds > 1 & inds < length(X));
        inds = inds(inds > 1 & inds < length(X));
    else
        % If there are no peaks, return an empty matrix
        freqs = [];
        inds = [];
    end
