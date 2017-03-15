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


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Process audio...
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sum to mono if necessary
    if size(x, 2) ~= 1
        x = sum(x,size(x,2))*0.5;
    end
    x = x';

    loopPoints(x, p);
end

