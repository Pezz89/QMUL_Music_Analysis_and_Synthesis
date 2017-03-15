function a = loopPoints(signal, p)
    % Calculate fundamental frequencies above a set inharmonicity threshold
    f0 = calcF0(signal, p);
    % Calculate energy of frames
    rms = RMS(signal, p);
    % Calculate transience of frames
    sf = spectralFlux(signal, p);

    s.sf = sf;
    s.rms = rms';
    s.f0 = f0;
    save('./analysis.mat', '-struct', 's');

    % Find frames where the standard deviation of f0, spectral flux and energy is below a
    % certain threshold.
    minNoPeriods = 20;

    % Consider only frames that are less than 25dB less than the maximum RMS
    % frame
    silenceThresh = 10^-2.5*max(rms);
    silenceMask = rms > silenceThresh;
    inharmonicMask = f0 == 0;

    % Calculate the standard deviation of consecutive frames to measure the
    % spread of values over time
    sfstd = movstd(sf, 3);
    rmsstd = movstd(rms, 3);
    f0std = movstd(f0, 3);
    % Inharmonic frames are set to the maximum standard deviation +10% as they
    % are most likely to be worse for looping than harmonic frames
    f0std(inharmonicMask) = max(f0std) + max(f0std)*0.1;

    % Set threshold initially at mix(feature) + 0.05 x min(feature) standard
    % deviation
    sfInc = 0.05 * min(sfstd(silenceMask));
    rmsInc = 0.05 * min(rmsstd(silenceMask));
    f0Inc = 0.05 * min(f0std(silenceMask));
    f0thresh = min(f0std(silenceMask));
    rmsthresh = min(rmsstd(silenceMask));
    sfthresh = min(sfstd(silenceMask));
    loopFound = false;

    while ~loopFound
        f0thresh = f0thresh + f0Inc;
        rmsthresh = rmsthresh + rmsInc;
        sfthresh = sfthresh + sfInc;
        % Find frames above thresholds for each feature
        f0candidates = find(f0std < f0thresh);
        rmscandidates = find(rmsstd < rmsthresh);
        sfcandidates = find(sfstd < sfthresh);

        f0CDiff = diff(f0candidates) == 1
        if isempty(f0CDiff)
            continue;
        end

        % Code adapted from https://uk.mathworks.com/matlabcentral/newsreader/view_thread/151318
        krn=[1 -1];
        changes=conv(krn, f0CDiff);

        % Calculate start and end window indexes of transient segments
        t_s = f0candidates(changes==1);
        t_e = f0candidates(changes==-1);

        f0MeanSeg = [];
        % Get mean f0 for each segment
        for i=1:length(t_s)
            f0MeanSeg = [f0MeanSeg, mean(f0(t_s(i):t_e(i)))];
        end
        % if segment is longer than the minimum number of periods at the
        % segment's mean f0 size, return it's start point's next zero crossing
        % and the integer multiple of the period closest to it's endpoint


        if any(f0MeanSeg)
            keyboard
        end

        % Convert window index to samples
        %transience_s = t_s * hopSize;
        %transience_e = t_e * hopSize;
    end


    % Calculate start and end times for segments of consecutive frames
    % Calculate the lengths of the start and end times
    % Find the largest segment
    % If no segment can be found, raise the threshold incremently until a
    % suitably long segment is found.

