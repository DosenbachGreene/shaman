function x = gencorr(y, r, x0)
% Generate a variable x whose correlation with y is r.  Takes an optional
% vector x0 with the same dimensions as y.  If x0 is not provided then
% defaults to a vector drawn from the standard normal distribution.

assert(isvector(y));
y = y(:);
assert(isscalar(r) && r >= -1 && r <= 1);
if nargin == 2
    x0 = normrnd(0, 1, length(y), 1);
end

resid = x0 - (y\x0)*y; % x0 - pinv(y)*x0*y;
x = r * std(resid) * y + resid .* std(y) * sqrt(1 - r^2);

end