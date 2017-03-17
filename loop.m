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
    y = zeros(1, length(signal) + (loopCount-1)*loopLength);

    y(1:loopStart-1) = signal(1:loopStart-1);

    if p.pyplot
        figure
        plot(signal(loopStart-40:loopEnd+40))
        hold on
        plot([zeros(1, 40), signal(loopStart:loopEnd-1), zeros(1, 41)])
    end

    loop = signal(loopStart:loopEnd-1);
    for i=0:loopCount-1
        y(loopStart+(loopLength*i):loopEnd-1+(loopLength*i))=y(loopStart+(loopLength*i):loopEnd-1+(loopLength*i))+loop;
    end
    y(loopEnd+(loopLength*i):end) = signal(loopEnd:end);

