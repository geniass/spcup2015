function [ y_new ] = preprocess( y, Fs, Fs_new )
%PREPROCESS Summary of this function goes here
%   Detailed explanation goes here

    y_new=y(1:min(length(y),100*Fs));
    % decimate
    Fs_old = Fs;
    Fs = Fs_new;
    y_new = resample(y_new, Fs, Fs_old);

    % bandpass filter
    b = fir1(50,[40,70].*2/Fs);
    y_new = filtfilt(b,1,y_new);

end

