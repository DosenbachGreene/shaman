function m = corrmat_unvectorize(v)
% Reconstitute a correlation matrix that was vectorized by corrmat_vectorize().
% 
% Example:
%     v = [5;7;8];
%     m = corrmat_unvectorize(v)
%     m =
%     
%         0     5     7
%         5     0     8
%         7     8     0

% Calculate matrix size
msize = (1 + sqrt(1 + 8*length(v))) / 2;
assert(mod(msize,1) == 0);
assert(msize > 1);

% Allocate matrix
m = zeros(msize, msize);

% Reconstitute lower triangular matrix
mcol = 1;
vrow = 1;
while mcol < msize
    i = mcol + 1;
    j = vrow + msize - mcol - 1;
    m(i:msize, mcol) = v(vrow:j);
    mcol = i;
    vrow = j + 1;
end

% Mirror lower triangular matrix to upper
m = tril(m)' + m;
