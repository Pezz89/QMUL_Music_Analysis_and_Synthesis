function analysis = spectralFlux(in, p)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Function to calculate Spectral Flux analysis for an input.
    % Inputs:
    %   in: Input audio vector
    %   wsize: Analysis window size
    %   hop: Analysis hop size
    % Returns:
    %   analysis: Normalised Spectral Flux analysis frames
    %   winCount: The total number of windows used in analysis of the input
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~isfield(p, 'wsize'); error('Window size not provided to spectral flux function'); end
    if ~isfield(p, 'hop'); p.hop=p.wsize; end

    wsize = p.wsize;
    hop = p.hop;
    % Allocate memory to store the current grain to be analysed
    grain = zeros(wsize,1);
    % Allocate memory to store the previous window's magnitude during analysis
    w1 = hanning(wsize); % Analysis Hanning window of length wsize
    mag1 = zeros(wsize/2,1);
    % Calculate the total number of windows to be analysed
    winCount = floor((length(in)-wsize)/hop);
    % Allocate memory to store outut analysis values
    analysis = zeros(winCount, 1);

    % Declare start and end indexes for reading from and writing to memory
    pin  = 0;
    pout = 1;
    pend = length(in)-wsize;
    % For each window in the source audio...
    while pout<winCount
        grain = in(pin+1:pin+wsize).* w1;
        f = fft(grain);
        % Calculate the magnitude of all non-mirrored FFT bins
        mag = abs(f(1:wsize/2))';

        % Calculate spectral flux analysis for the current window
        analysis(pout) = sqrt(sum((mag-mag1).^2))/(wsize/2);

        %%%%% Alternate method for calculating Spectral Flux... %%%%%
        %
        % mag_diff = mag-mag1
        % analysis(pout) = sum(mag_diff-abs(mag_diff)/2);
        %
        %%%%%

        % Store magnitude of current frame for use in the next frame
        mag1 = mag;

        % Increment read and write pointers
        pin  = pin + hop;
        pout = pout + 1;
    end

