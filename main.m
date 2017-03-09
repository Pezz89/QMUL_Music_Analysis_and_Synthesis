function main()
    midiFileName = './media/bach_pr1.mid';
    wav1FileName = './media/bach_pr1_a.wav';
    wav2FileName = './media/bach_pr1_b.wav';
    f = load('./media/notes.mat');
    transcriptionMatrix = f.notes

    [in, fs] = wavread(audioFile);
    [FFTidx, Fp_est, Fp_corr] = find_pitch_fft(in, blackmanharris(4096), 4096, fs, 4096/2, 20.0, 20000.0, -20.0)

    % [f0 B] = fundamentals(wav1FileName, transcriptionMatrix)
