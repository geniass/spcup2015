function [ F ] = freq_est_music( y, Fs )
%FREQ_EST estimates the single instantaneous frequency in the vector y
% returns the vector of instantaneous frequencies



W = 1*Fs;
L = 1;

N = length(y);
F = zeros(1,length(1:(N-W)));

m = window_mat(y, L, W)';
M_1 = size(m,1);

clear y

window = hamming(W);
% padding = zeros(b*W,1);

F = zeros(1,length(1:(N-W)));

tic
% matlabpool(3)

%  profile('on', '-history', '-detail', 'builtin','-memory')
parfor n = 1:size(m,2)
    [~,R] = corrmtx(m(:,n),24,'modified');
    [W,~] = rootmusic(R,2,Fs,'corr');
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

% figure
% t = 0:1/Fs:(length(F)-1)/Fs;
%  plot(t,F)
end