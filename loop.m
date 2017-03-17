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
    % Allocate memory for output
    y = zeros(1, length(signal) + (loopCount-1)*loopLength);

    % Copy all audio before the loop point
    y(1:loopStart-1) = signal(1:loopStart-1);

    % Get loop frames between start and end segments
    loop = signal(loopStart:loopEnd-1);
    % Repeatedly copy loop to output at incrementing frames based on the number
    % of loops requested
    for i=0:loopCount-1
        y(loopStart+(loopLength*i):loopEnd-1+(loopLength*i))=y(loopStart+(loopLength*i):loopEnd-1+(loopLength*i))+loop;
    end
    % Copy all other audio to end of buffer
    y(loopEnd+(loopLength*i):end) = signal(loopEnd:end);

    % Plot zero crossing example data using python
    if p.pyplot
        n = 300
        s.raw = signal(loopStart-n:loopEnd+n);
        s.clipped = [zeros(1, n), signal(loopStart:loopEnd-1), zeros(1, n+1)];
        save('./analysis3.mat', '-struct', 's');
    end
