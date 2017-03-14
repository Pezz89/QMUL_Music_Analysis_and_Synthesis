function main()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % User settings
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    audioFile = './media/vibrato.wav';
    %Load audiofile
    [x, FS] = audioread(audioFile);

    % Set samplerate
    p.FS = FS;
    % set f0 minimum and maximum period thresholds
    p.f0min=50;
    % Set f0 and RMS window and hop sizes
    p.wsize = 2048;
    p.hop = round(p.wsize);

    % Set F0 harmonicity threshold
    p.fDelta = 0;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Process audio...
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sum to mono
    %x = sum(x,size(x,2))*0.5;
    x = x';

    loopPoints(x, p);
end

