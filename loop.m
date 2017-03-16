function y = loop(signal, duration, loopPoints, p)
    % Get loop start and end points
    loopStart = loopPoints(1);
    loopEnd = loopPoints(2);
    loopLength = loopEnd - loopStart;

    if duration < signal
        error('Duration is shorter than the original signal')
    end
    %
    % Calculate difference in samples between the current duration and the new
    % duration
    durationDiff = duration - length(signal);
    % Calculate the number of times that the loop will fully fit into this
    loopCount = floor(durationDiff / loopLength);

    y = signal(1:loopStart);

    figure
    plot(signal(loopStart-40:loopEnd+40))
    hold on
    plot([zeros(1, 40), signal(loopStart:loopEnd-1), zeros(1, 41)])

    for i=1:loopCount
        y = [y, signal(loopStart:loopEnd-1)];
    end
    y = [y, signal(loopEnd+1:end)];

