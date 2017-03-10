function [freqs, inds] = findPeaks(X, frequencyPerBin)
    Xmax = max(X);
    ThreshInds = find(X > (movmean(X, 200) + (0.5 * movstd(X, 200))) & X > (10^-2.5)*Xmax);

    if false
        figure
        plot((movmean(X, 200) + (0.5 * movstd(X, 200))), 'g')
        hold on
        plot(X)
        hold on
    end
    if any(ThreshInds)
        peaks = X(ThreshInds);
        [peaks, locs] = findpeaks(X);
        [~, ~, ib] = intersect(ThreshInds, locs);
        peaks = peaks(ib);
        locs = locs(ib);
        freqs = locs * frequencyPerBin;

        ampfreq = [peaks freqs locs];

        sorted = sortrows(ampfreq,1);

        freqs = flip(sorted(:,2));
        inds = flip(sorted(:,3));
        freqs = freqs(inds > 1 & inds < length(X));
        inds = inds(inds > 1 & inds < length(X));
    else
        freqs = [];
        inds = [];
    end
