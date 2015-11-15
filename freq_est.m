clear all

dir = '/home/ari/development/spcup/Grid_B/Power_recordings/';
file = 'Train_Grid_B_P1.wav';
[y,Fs] = wavread(strcat(dir,file));
y = y(1:1000*Fs);

%  [y,n,window_freqs,sample_freqs]=synth_enf_signal(Fs, n_windows, window_size);
% n = 0:1/Fs:4;
% y = cos(2*pi*50*n.*(n<2) + 2*pi*55.1*n.*(n>=2));
% plot(n,y)
size(y)

% decimate
Fs_old = Fs;
Fs = 300;
y = resample(y, Fs, Fs_old);

b = fir1(50,[40,70].*2/Fs);
y = filtfilt(b,1,y);

W = 1*Fs;
L = 1;

N = length(y)
F = zeros(1,length(1:(N-W)));

m = window_mat(y, L, W)';
M_1 = size(m,1)

clear y

window = hamming(W);
% padding = zeros(b*W,1);

F = zeros(1,length(1:(N-W)));

tic
% matlabpool(3)

%  profile('on', '-history', '-detail', 'builtin','-memory')
parfor n = 1:size(m,2)
    [~,R] = corrmtx(m(:,n),24,'modified');
    [W,P] = rootmusic(R,2,Fs,'corr');
    F(n) = W(1);
end
toc

% N = length(y)
% wind_half_size = window_size/4;
% F = zeros(1,N);
% tic
% parfor n = 1:N
%     [X,R] = corrmtx(y(max(1,n-wind_half_size) : min(N, n + wind_half_size)),24,'modified'); % Why 12? ... ?
%     [W,P] = rootmusic(X,2,Fs);
%     F(n) = W(1);
% end

% toc

figure
t = 0:1/Fs:(length(F)-1)/Fs;
 plot(t,F)
%plot(t,sample_freqs,t,F)

% runs out of memory
% figure
% level = 7;
% wpt = wpdec(y,level,'sym8');
% [Spec,Time,Freq] = wpspectrum(wpt,Fs,'plot');


sprintf('Mean: %f\nVar: %f\n',mean(F),var(F));
