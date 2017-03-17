function [f0, f0Trans] = calcF0(x, p)
    if ~isfield(p, 'fDelta'); error('fDelta not provided to calcF0'); end
    % Calculate f0 of signal frames filtered by inharmonicity
    f = yin(x, p);
    f0 = f;

    % Apply further processing only to values that are harmonic. Avoids divide
    % by zero error
    f0Trans = zeros(size(f));
    mask = f > 0;

    % Map to semi-tone pitch spectrum
    f(mask) = 12*log2(f(mask)/440.0)+69.0;

    % place data roughly around the center of the chroma spectrum
    f = f - (median(f(mask))+6);

    % Wrap octaves to remove octave error (and normalise to range of 0-1)
    f0Trans(mask) = mod(f(mask), 12)/12;
