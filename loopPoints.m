function a = loopPoints(signal, p)
    % Calculate fundamental frequencies above a set inharmonicity threshold
    f0 = calcF0(signal, p);
    % Calculate energy of frames
    rms = RMS(signal, p);
    % Calculate transience of frames
    sf = spectralFlux(signal, p);

    s.sf = sf;
    s.rms = rms';
    s.f0 = f0;
    save('./analysis.mat', '-struct', 's');

    % Find frames where the standard deviation of f0, transience and energy is below a
    % certain threshold.
    % If no segment can be found, lower the threshold incremently until a
    % suitably long segment is found.

