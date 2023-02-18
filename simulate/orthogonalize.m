function A = orthogonalize(A,B)
% make A orthogonal to B

X = [ones(size(A,1),1), B];
A = A - X*(X\A);

end