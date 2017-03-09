function main2()
    clear all;
    close all;

    bach_pr1;

    [Y, FS] = audioread('./media/bach_pr1_a.wav');
    Ymono = sum(Y,size(Y,2));
    %%
    figure;
    plot(Ymono);
    hold on;
    stem(notes(:,1)*FS,notes(:,3).*1/128, 'r');
    axis([5*FS 20*FS -1 1]);
    %%
    Wlen = 4096;
    w = blackmanharris(Wlen);

    sample =10*FS;
    grain = Ymono(sample+1:sample+Wlen) .* w;
    plot(grain);

    %%
    f     = fft(fftshift(grain));
    X     = abs(f);
    phi   = angle(f);

    figure
    peaks = X > mvav(X, 200)' + (0.5 * mvstd(X, 200)');
    stem(peaks, 'r')
    hold on
    plot(X/max(X))
    %%
    %plot(X([1:Wlen/40]));
    % hold on;
    % stem(X([1:Wlen/40]));
    % %%
    % [p, l] = findpeaks(X([1:Wlen/2]));
    % hold on;
    % stem(X([1:Wlen]));
    % %i = find(p==max(p));
    % i2 = find(X==max(X),1);
    % %stem(l(i2), max(p));
    % stem(i2, max(X), 'g');
    % stem(i2-1, X(i2-1), 'g');
    % stem(i2+1, X(i2+1), 'g');
    % axis([1 Wlen/40 0 60]);
    % %%
    % alpha = X(i2-1);
    % beta = X(i2);
    % gamma = X(i2+1);
    %
    % p = 0.5*((alpha-gamma)/(alpha-2*beta+gamma));
    % peak = i2+p;
    % %%
    % peak_magnitude = beta-0.25*(alpha-gamma)*p;
    % %%
    % hold on;
    % stem(peak, peak_magnitude, 'r');
end

function output = mvav(x,n)
% Preallocate output
output=NaN(1,numel(x));
% Find mid point of n
L = length(x);
midPoint = round(n/2);
x = [zeros(midPoint, 1); x; zeros(midPoint, 1)];
% For length of input, minus n
    for a = midPoint:L

        % Find index range to take average over
        b=a+midPoint;
        % Calculate mean
        output(a-midPoint+1) = mean(x(a-midPoint+1:b));

    end
end

function output = mvstd(x,n)
% Preallocate output
output=NaN(1,numel(x));
% Find mid point of n
L = length(x);
midPoint = round(n/2);
x = [zeros(midPoint, 1); x; zeros(midPoint, 1)];
% For length of input, minus n
    for a = midPoint:L

        % Find index range to take average over
        b=a+midPoint;
        % Calculate mean
        output(a-midPoint+1) = std(x(a-midPoint+1:b));

    end
end
