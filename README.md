# spcup2015

First thing to do is extract the ENF signal from recordings:

1) freq_est.m goes through every file in the training dataset and calls preprocess.m. This downsamples the data (much faster) and selects only a small part of the data (otherwise the script will takes days to run). It then bandpass filters the data
2) calls freq_est_stft or freq_est_music to extract the ENF signal
   freq_est_music uses the rootMUSIC algorithm (don't really understand how it works. But it does. And it's really slow)
