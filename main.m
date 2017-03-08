function main()
    midiFileName = './media/bach_pr1.mid';
    wav1FileName = './media/bach_pr1_a.wav';
    wav2FileName = './media/bach_pr1_b.wav';
    f = load('./media/notes.mat');
    transcriptionMatrix = f.notes

    [f0 B] = fundamentals(wav1FileName, transcriptionMatrix)
