function loopInds = loopPoints(signal, p)
    % Calculate fundamental frequencies above a set inharmonicity threshold
    [f0, f0Trans] = calcF0(signal, p);
    % Calculate energy of frames
    rms = RMS(signal, p);
    % Calculate transience of frames
    sf = spectralFlux(signal, p);

    if p.pyplot
        s.sf = sf;
        s.rms = rms';
        s.f0 = f0;
        save('./analysis.mat', '-struct', 's');
        [status,cmdout] = system('./GeneratePlots.py');
    end

    if ~isfield(p, 'minPeriod'); p.minPeriod=350; end
    % Find frames where the standard deviation of f0, spectral flux and energy is below a
    % certain threshold.
    minNoPeriods = p.minPeriod;

    % Consider only frames that are less than 25dB less than the maximum RMS
    % frame
    silenceThresh = 10^-2.5*max(rms);
    %silenceMask = rms > silenceThresh;
    silenceMask = rms > silenceThresh;
    inharmonicMask = f0Trans == 0;

    % Calculate the standard deviation of consecutive frames to measure the
    % spread of values over time
    sfstd = movstd(sf, 9);
    rmsstd = movstd(rms, 9);
    f0Transstd = movstd(f0Trans, 9);
    if p.pyplot
        s.sf = sfstd;
        s.rms = rmsstd';
        s.f0 = f0Transstd;
        save('./analysis.mat', '-struct', 's');
        [status,cmdout] = system('./GeneratePlots.py');
    end
    % Inharmonic frames are set to the maximum standard deviation +10% as they
    % are most likely to be worse for looping than harmonic frames
    f0Transstd(inharmonicMask) = max(f0Transstd) + max(f0Transstd)*0.1;

    p.thresholdInc = 0.0001;
    sfInc = p.thresholdInc * (max(sfstd(silenceMask)) - min(sfstd(silenceMask)));
    rmsInc = p.thresholdInc * (max(rmsstd(silenceMask)) - min(rmsstd(silenceMask)));
    f0TransInc = p.thresholdInc * (max(f0Transstd(silenceMask)) - min(f0Transstd(silenceMask)));
    f0Transthresh = min(f0Transstd(silenceMask));
    rmsthresh = min(rmsstd(silenceMask));
    sfthresh = min(sfstd(silenceMask));
    loopFound = false;

    while ~loopFound
        f0Transthresh = f0Transthresh + f0TransInc;
        rmsthresh = rmsthresh + rmsInc;
        sfthresh = sfthresh + sfInc;
        if f0Transthresh > max(f0Transstd)*1.4 || rmsthresh > max(rmsstd)*1.4 || ...
            sfthresh > max(sfstd)*1.4
            error('No segments found, try lowering the minimum period count or increment rate');
        end
        % Find frames below thresholds for each feature
        f0TransCandidates = find(f0Transstd < f0Transthresh);
        rmscandidates = find(rmsstd < rmsthresh);
        sfcandidates = find(sfstd < sfthresh);

        % Find consecutive frames
        f0TransCDiff = diff(f0TransCandidates) == 1;
        if isempty(f0TransCDiff)
            continue;
        end

        % Code adapted from https://uk.mathworks.com/matlabcentral/newsreader/view_thread/151318
        krn=[1 -1];
        changes=conv(krn, f0TransCDiff);

        % Calculate start and end window indexes of transient segments
        t_s = f0TransCandidates(changes==1);
        t_e = f0TransCandidates(changes==-1);

        f0AvrSeg = [];
        % Get mean f0Trans for each segment
        for i=1:length(t_s)
            f0AvrSeg = [f0AvrSeg, median(f0(t_s(i):t_e(i)))];
        end
        % if segment is longer than the minimum number of periods at the
        % segment's mean f0 size, return it's start point's next zero crossing
        % and the integer multiple of the period closest to it's endpoint

        % Convert window index to samples
        segStart = (t_s * p.hop) + round(p.wsize/2);
        segEnd = (t_e * p.hop) + round(p.wsize/2);
        segLength = segEnd - segStart;
        viableSegs = (((1./f0AvrSeg)*p.FS)*minNoPeriods) < segLength;
        if any(viableSegs)
            segLength(~viableSegs) = Inf;
            [~, ind] = min(segLength);
            finalStart = segStart(ind);
            finalEnd = segEnd(ind);
            finalMF0 = f0AvrSeg(ind);
            finalLength = segLength(ind);
            loopFound = true;
        end
    end

    if finalLength == 0
        error('Segment length is 0. This shouldnt happen...')
    end

    % Refine start and end points based on zero-crossings
    % Find the nearest zero crossing to the middle of the start point's window
    % Code adapted from: https://uk.mathworks.com/matlabcentral/newsreader/view_thread/48430
    halfWindow = round(p.wsize/2);
    x = signal(finalStart-halfWindow:finalStart+halfWindow-1);
    signum = sign(x);    % get sign of data
    signum(x==0) = 1;   % set sign of exact data zeros to positiv
    zeroX = find(diff(signum)~=0);  % get zero crossings by diff ~= 0

    tmp = abs(zeroX-halfWindow);
    [~, idx] = min(tmp); %index of closest value
    finalStart = finalStart-halfWindow+zeroX(idx);
    % Find the nearest sample index to the end index that is an integer multiple of the mean f0
    % period
    % Calculate the final period in samples
    finalPeriod = (1./finalMF0)*p.FS;

    % Calculate index of the nearest multiple to the fundamental period
    finalEnd = floor(finalLength/finalPeriod)*(finalPeriod)+finalStart;
    % Find the nearest zero crossing to the index found
    x = signal(finalEnd-halfWindow:finalEnd+halfWindow-1);
    signum = sign(x);    % get sign of data
    signum(x==0) = 1;   % set sign of exact data zeros to positiv
    zeroX = find(diff(signum)~=0);  % get zero crossings by diff ~= 0

    tmp = abs(zeroX-halfWindow);
    [~, idx] = min(tmp); %index of closest value
    finalEnd = finalEnd-halfWindow+zeroX(idx);

    if finalEnd - finalStart == 0
        error('Segment length is 0. This shouldnt happen...')
    end
    loopInds = [round(finalStart), round(finalEnd)];

    % Calculate start and end times for segments of consecutive frames
    % Calculate the lengths of the start and end times
    % Find the largest segment
    % If no segment can be found, raise the threshold incremently until a
    % suitably long segment is found.

