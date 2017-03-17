function main()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % User settings
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    audioFile = './media/sguitar.aiff';
    %Load audiofile
    [x, FS] = audioread(audioFile);

    p.ratio = 1.5;
    % Set samplerate
    p.FS = FS;
    % set f0 minimum and maximum period thresholds
    p.f0min=50;
    % Set f0 and RMS window and hop sizes
    p.wsize = 2048;
    p.hop = round(p.wsize*0.5);

    % Set F0 harmonicity threshold
    p.fDelta = 0;

    % Set the minimum number of f0 periods that is allowed in a segment
    p.minPeriod = 14;

    % Plot results using python script
    p.pyplot = true;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Process audio...
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sum to mono if necessary
    if size(x, 2) ~= 1
        x = sum(x,size(x,2))*0.5;
    end
    % Normalise audio
    x = x/max(x);
    x = x';

    % Calculate start and end indexes in samples for loop
    loopInds = loopPoints(x, p);

    % Plot start and end indexes using python
    if p.pyplot
        s.loopInds = loopInds;
        save('./loopInds.mat', '-struct', 's');
        [status,cmdout] = system('./GeneratePlots.py');
    end

    % Apply looping to stretch audio to required duration
    y = loop(x, length(x) * p.ratio, loopInds, p);
    audiowrite('./out.wav', y, p.FS);
end

