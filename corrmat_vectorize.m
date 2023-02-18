function v = corrmat_vectorize(m)
% CORRMAT_VECTORIZE  Flatten a correlation matrix into a column vector of
% its lower triangular matrix.
% 
% Example:
%     m = [1,2,3;4,5,6;7,8,9]
%     m =
%     
%         1     2     3
%         4     5     6
%         7     8     9
%     
%     v = corrmat_vectorize(m)
%     v =
%     
%         4     7     8
% 
% See also: CORRMAT_UNVECTORIZE, TRIL

% Matrix must be square.
msize = size(m);
assert(msize(1) == msize(2));

% Allocate output vector v.
msize = msize(1);
v = zeros(1, (msize*msize-msize)/2);

% Extract lower triangular matrix and vectorize into v.
% Each column in the lower triangular matrix becomes the next segment in v.
mcol = 1;
vcol = 1;
while mcol < msize
    i = mcol + 1;
    j = vcol + msize - mcol - 1;
    v(vcol:j) = m(i:msize, mcol);
    mcol = i;
    vcol = j + 1;
end
