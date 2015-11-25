clear all

matlab_dir = pwd;
data_dir = '/home/ari/development/spcup/';
output_dir = '/home/ari/development/spcup/training_data';
audio_dir = 'Audio_recordings';
power_dir = 'Power_recordings';
Fs_new = 300;
cd(data_dir);
grids = dir('Grid*');


% run frequency estimation on every file
for grid = grids'        % must have ' for some reason ?!?
    grid_dir = fullfile(data_dir, grid.name);
    
    % Audio
    cd(fullfile(grid_dir,audio_dir));
    files = dir('*.wav');
    for file = files'
        fprintf('%s\n', file.name)
        % process
        filename = fullfile(data_dir,grid.name,audio_dir,file.name);
%         filename = '/home/ari/development/spcup/Grid_G/Power_recordings/Train_Grid_G_P1.wav';
        [y,Fs] = wavread(filename);
        
        cd(matlab_dir);
        y = preprocess(y, Fs, Fs_new);
        F = freq_est_music(y, Fs_new);
        save(fullfile(output_dir,strcat(file.name,'.mat')),'F');
        
        clear y
        clear F
    end
    
    % Power
    cd(fullfile(grid_dir,power_dir));
    files = dir('*.wav');
    for file = files'
        fprintf('%s\n', file.name)
        % process
        filename = fullfile(data_dir,grid.name,power_dir,file.name);
        filename = '/home/ari/development/spcup/Grid_G/Power_recordings/Train_Grid_G_P1.wav';
        [y,Fs] = wavread(filename);
        
        cd(matlab_dir);
        y = preprocess(y, Fs, Fs_new);
        F = freq_est_music(y, Fs_new);
        save(fullfile(output_dir,strcat(file.name,'.mat')),'F');
        
        clear y
        clear F
    end
end



fprintf('DONE!\n')
pause



file = 'Train_Grid_B_P1.wav';
[y,Fs] = wavread(strcat(data_dir,file));
 y = y(1:1000*Fs);
size(y)
%  [y,n,window_freqs,sample_freqs]=synth_enf_signal(Fs, n_windows, window_size);
% n = 0:1/Fs:4;
% y = cos(2*pi*50*n.*(n<2) + 2*pi*55.1*n.*(n>=2));
% plot(n,y)







% decimate
Fs_old = Fs;
Fs = 300;
y = resample(y, Fs, Fs_old);
size(y)

% bandpass filter
b = fir1(50,[40,70].*2/Fs);
y = filtfilt(b,1,y);

window_size = 1000;
M = floor(window_size / 3);
n_windows = 20;

% frame length = D x L
L = 1;   % hop size (0.1s)
% b = 8;      % zero-padding factor
D = 1*Fs;      % sample interval mult. factor


N = length(y)
F = zeros(1,length(1:(N-D*L)));

% window and zero-padding for better FFT resolution
P = 2^16;       % FFT size (power of 2)
window = blackman(D*L);
% padding = zeros(8192 - D*L,1);
m = window_mat(y, L, D*L)';
m = [m .* repmat(window,1, size(m,2))];% ; repmat(padding,1,size(m,2))];

clear y
clear window
clear padding

fprintf('Starting loop...\n')
tic
% matlabpool(3)

% profile off
%  profile('on', '-history', '-detail', 'builtin','-memory')
parfor n = 1:size(m,2)
%     x = [m(:,n) .* window ; padding];
%     x = [x ; padding];   % zero-padding
    x = m(:,n);
    
    Y = fft(x,P);
    P2 = abs(10*Y/D/L);
    P1 = P2(1:P/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    f = Fs*(0:(P/2))/P;
    % https://ccrma.stanford.edu/~jos/parshl/Peak_Detection_Steps_3.html#sec:peakdet
    [~,i] = max(P1);                    % find max power at i
    b = 20*log10(P1(i));
    c = 20*log10(P1(i+1));
    a = 20*log10(P1(i-1));             % find two adjacent bin magnitudes
    
    p = 0.5*(a-c)/(a-2*b+c);   
    bin = round(i + p);
    f_max = Fs*bin/P;
    F(n) = Fs*i/P;
    
%     if mod(n,100)==0
%         sprintf('Completion: %i %%\n', round(n*100/(N-D*L)))
%         sprintf('vals: %f, %f, %f\n',P1(i-1),P1(i),P1(i+1))
%         sprintf('max: %f, interp: %f\n', i, bin)
%     end
end
%  profile viewer
% profile stop
toc



figure
t = 0:L/Fs:(length(F)-1)/Fs;
 plot(t,F)
%plot(t,sample_freqs,t,F)

% runs out of memory
% figure
% level = 7;
% wpt = wpdec(y,level,'sym8');
% [Spec,Time,Freq] = wpspectrum(wpt,Fs,'plot');


sprintf('Mean: %f\nVar: %f\n',mean(F),var(F));
