function pitches = removeHarmonics(freq, relativeError)

pitches = [];

if ~isempty(freq)   
    freqs = [ 1:1:length(freq) ; freq ];
    freqs = sortrows(freqs',2)';
    while ~isempty(freqs) 
       currentf  = freqs(2,1);
       pitches = [pitches freqs(:,1);]; 
    %    remainder = mod((freqs / currentf),1);
       remainder = mod(freqs(2,:), currentf);
       allowedError = relativeError * freqs(2,:);
       hid = find( remainder > allowedError & abs(remainder-currentf) > allowedError);
       freqs = freqs(:,hid');
%        freqs = freqs(:, find( currentf ~= freqs(2,:) ) );
    end

    sortrows(pitches,1);
    pitches = pitches(2,:);
end