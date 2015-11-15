function [ m ] = window_mat( data, hop, width )
%WINDOW_MAT Builds a matrix where each row is the output of a sliding
%window applied to data
% NOTE IT IS CURRENTLY BROKEN FOR HOPS OF MORE THAN 2

N = length(data);
M = floor((N-width)/hop + 1);
m = zeros(M, width);
for n = 1:M
    m(n,:) = data(n*hop - floor(hop/2) : n*hop - floor(hop/2) + width - 1);
end

end

