function r = RMS(signal, p)
    if ~isfield(p, 'maxprd'); p.maxprd=100; end
    if ~isfield(p, 'wsize'); p.wsize=p.maxprd; end
    if ~isfield(p, 'hop'); p.hop=p.wsize; end
    if ~isfield(p, 'FS'); error('Samplerate not provided to RMS function'); end

    if min(size(signal)) ~= 1; error('data should be 1D'); end
    x=signal(:)/max(signal(:));
    nsamples=numel(x);

    nframes=floor((nsamples-p.wsize)/p.hop);
    r=zeros(1,nframes);

    for k=1:nframes
        % Get current frame of audio
        start=(k-1)*p.hop; % offset of frame
        xx=x(start+1:start+p.wsize,:);

        % Turns out theres a built in matlab function for this now...
        % would have been: r(k) = sqrt(mean(xx.^2));
        r(k) = rms(xx);
    end
