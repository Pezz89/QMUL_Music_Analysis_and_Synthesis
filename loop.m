function y = loop(signal, duration)
    % Set a size in samples for the cross fade between consecutive loops
    crossfadeSize = 100*FS
    % Last loop must be truncated from thebegining so that the end matches up
    % with the non loop bit...
