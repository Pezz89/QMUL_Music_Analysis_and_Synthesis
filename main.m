function main()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % User settings
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    audioFile = './media/sguitar.aiff';
    %Load audiofile
    [x, FS] = audioread(audioFile);

    % Set samplerate
    p.FS = FS;
    % set f0 minimum and maximum period thresholds
    p.f0min=50;
    % Set f0 and RMS window and hop sizes
    p.wsize = 2048;
    p.hop = round(p.wsize*0.5);

    % Set F0 harmonicity threshold
    p.fDelta = 0;

    p.minPeriod = 100;

    % Plot results using python script
    p.pyplot = false;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Process audio...
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sum to mono if necessary
    if size(x, 2) ~= 1
        x = sum(x,size(x,2))*0.5;
    end
    x = x';

    % Calculate start and end indexes in samples for loop
    loopInds = loopPoints(x, p);

    if p.pyplot
        s.loopInds = loopInds;
        save('./loopInds.mat', '-struct', 's');
    end

    y = loop(x, length(x) * 1.5, loopInds);
    audiowrite('./out.wav', y, p.FS);
end

