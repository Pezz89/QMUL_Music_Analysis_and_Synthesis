function y = loop(signal, duration, loopPoints, p)
    % Set crossfade size
    crossFade = round(0.05*p.FS)
    % Get loop start and end points
    loopStart = loopPoints(1);
    loopEnd = loopPoints(2);
    loopLength = loopEnd - loopStart;
    if crossFade > loopLength / 2
        crossFade = floor(loopLength / 2);
    end

    if duration < signal
        error('Duration is shorter than the original signal')
    end
    %
    % Calculate difference in samples between the current duration and the new
    % duration
    durationDiff = duration - length(signal);
    % Calculate the number of times that the loop will fully fit into this
    loopCount = floor(durationDiff / (loopLength-crossFade*2));
    loopIncrement = (loopLength-(crossFade*2));
    y = zeros(1, duration);

    fadeIn = sqrt(0.5*(1+linspace(-1, 1, crossFade)))
    fadeOut = sqrt(0.5*(1-linspace(-1, 1, crossFade)))
    y(1:loopStart+crossFade) = signal(1:loopStart+crossFade);
    y(loopStart:loopStart+crossFade-1) = y(loopStart:loopStart+crossFade-1) .* fadeOut;

    figure
    plot(signal(loopStart-40:loopEnd+40))
    hold on
    plot([zeros(1, 40), signal(loopStart:loopEnd-1), zeros(1, 41)])
    loop = signal(loopStart:loopEnd-1);
    loop(1:crossFade) = loop(1:crossFade) .* fadeIn;
    loop(end-crossFade+1:end) = loop(end-crossFade+1:end) .* fadeOut;
    for i=0:loopCount-1
        y(loopStart+(loopIncrement*i):loopEnd-1+(loopIncrement*i))=y(loopStart+(loopIncrement*i):loopEnd-1+(loopIncrement*i))+loop;
    end
    %y = [y, signal(loopEnd+1:end)];

