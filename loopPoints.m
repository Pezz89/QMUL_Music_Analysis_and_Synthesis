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
        save('./analysis1.mat', '-struct', 's');
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

    % Calculate the absolute standard deviation of consecutive frames to
    % measure the spread of values over time
    sfstd = abs(movstd(sf, 3));
    rmsstd = abs(movstd(rms, 3));
    f0TransStd = abs(movstd(f0Trans, 3));

    % Plot standard deviation of features in python
    if p.pyplot
        s.sf = sfstd;
        s.rms = rmsstd';
        s.f0 = f0TransStd;
        save('./analysis2.mat', '-struct', 's');
        [status,cmdout] = system('./GeneratePlots.py');
    end

    % Inharmonic frames are set to the maximum standard deviation +10% as they
    % are most likely to be worse for looping than harmonic frames
    f0TransStd(inharmonicMask) = max(f0TransStd) + max(f0TransStd)*0.1;

    % Set an increment factor as a percentage of the range for each feature
    p.thresholdInc = 0.00005;
    sfInc = p.thresholdInc * (max(sfstd(silenceMask)) - min(sfstd(silenceMask)));
    rmsInc = p.thresholdInc * (max(rmsstd(silenceMask)) - min(rmsstd(silenceMask)));
    f0TransInc = p.thresholdInc * (max(f0TransStd(silenceMask)) - min(f0TransStd(silenceMask)));
    % Set thresholds to the minimum value of each of the features
    f0Transthresh = min(f0TransStd(silenceMask));
    rmsthresh = min(rmsstd(silenceMask));
    sfthresh = min(sfstd(silenceMask));
    loopFound = false;

    % Search for groups of frames that are long enough to form a segment based
    % on incrementing thresholds
    while ~loopFound
        % Increment thresholds
        f0Transthresh = f0Transthresh + f0TransInc;
        rmsthresh = rmsthresh + rmsInc;
        sfthresh = sfthresh + sfInc;
        % If the threshold is larger than the maximum values of all features
        % and a frame has not been found, error.
        % TODO: Set this to a less arbritrary number than 1.4...
        if f0Transthresh > max(f0TransStd)*1.4 || rmsthresh > max(rmsstd)*1.4 || ...
            sfthresh > max(sfstd)*1.4
            error('No segments found, try lowering the minimum period count or increment rate');
        end
        % Find frames below thresholds for each feature
        f0TransCandidates = find(f0TransStd < f0Transthresh);
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
            % Find the segment with the lowest accumulated standard deviation
            % accross all features
            vt_s = t_s(viableSegs)';
            vt_e = t_e(viableSegs)';
            RMSStdAvrSeg = [];
            f0StdAvrSeg = [];
            sFStdAvrSeg = [];
            for i=1:length(vt_s)
                RMSStdAvrSeg = [RMSStdAvrSeg, mean(rmsstd(vt_s(i):vt_e(i)))];
                f0StdAvrSeg = [f0StdAvrSeg, mean(f0TransStd(vt_s(i):vt_e(i)))];
                sFStdAvrSeg = [sFStdAvrSeg, mean(sfstd(vt_s(i):vt_e(i)))];
            end
            match = zeros(1, length(segLength));
            match(~viableSegs) = Inf;
            match(viableSegs) = RMSStdAvrSeg + f0StdAvrSeg + sFStdAvrSeg;
            % Choose the segment with the lowest accumulated average standard
            % deviation
            [~, ind] = min(match);
            % Get data for chosen segment
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
    % Find the nearest positive zero crossing to the middle of the start point's window
    % Code adapted from: https://uk.mathworks.com/matlabcentral/newsreader/view_thread/48430
    halfWindow = round(p.wsize/2);
    x = signal(finalStart-halfWindow:finalStart+halfWindow-1);
    signum = sign(x);    % get sign of data
    signum(x==0) = 1;   % set sign of exact data zeros to positiv
    zeroX = find(diff(signum)>0);  % get zero crossings by diff ~= 0

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
    zeroX = find(diff(signum)>0);  % get zero crossings by diff ~= 0

    tmp = abs(zeroX-halfWindow);
    [~, idx] = min(tmp); %index of closest value
    finalEnd = finalEnd-halfWindow+zeroX(idx);

    if finalEnd - finalStart == 0
        error('Segment length is 0. This shouldnt happen...')
    end
    loopInds = [round(finalStart), round(finalEnd)];
