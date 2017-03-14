function f0 = calcF0(x, p)
    if ~isfield(p, 'fDelta'); error('fDelta not provided to calcF0'); end
    % Calculate f0 and inharmonicity of signal frames
    f = yin(x, p);

    f0 = zeros(size(f));
    mask = f > 0;
    f(mask) = 12*log2(f(mask)/440.0)+69.0;

    f = f - (median(f(mask))+6);

    f0(mask) = mod(f(mask), 12);
    % Calculate frequency in Hz from period
    %f0 = 1./(f.prd*p.FS);
    %f0 = f.prd;
    % Get inharmonicity factor for each frame
    %fHarm = f.ap;

    % Filter out any frequency estimates below the harmonicity threshold.
    %f0(fHarm < p.fDelta) = 0;
