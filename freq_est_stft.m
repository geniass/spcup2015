function [ F ] = freq_est_stft( y, Fs )
%FREQ_EST_STFT estimates the single instantaneous frequency in the vector y
% using the short time Fourier transform method (spectrogram)
% returns the vector of instantaneous frequencies


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
P = 2^17;       % FFT size (power of 2)
window = blackman(D*L);
% padding = zeros(8192 - D*L,1);
m = window_mat(y, L, D*L)';
m = [m .* repmat(window,1, size(m,2))];% ; repmat(padding,1,size(m,2))];

clear y
clear window
clear padding

fprintf('Starting loop...\n')
fprintf('No. windows: %i\n', size(m,2))
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
    
    % https://ccrma.stanford.edu/~jos/parshl/Peak_Detection_Steps_3.html#sec:peakdet
    [~,i] = max(P1);                    % find max power at i
    b = 20*log10(P1(i));
    c = 20*log10(P1(i+1));
    a = 20*log10(P1(i-1));             % find two adjacent bin magnitudes
    
    p = 0.5*(a-c)/(a-2*b+c);   
    bin = round(i + p);
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


end

