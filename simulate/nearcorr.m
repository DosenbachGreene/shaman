function X = nearcorr(A)
% Find the nearest correlation matrix to A.
% Based on the alternating projection algorithm by Nick Higham, see:
% https://nhigham.com/2013/02/13/the-nearest-correlation-matrix/
% And also:
%   N. J. Higham, 13/6/01, updated 30/1/13, 15/11/14, 07/06/15.
%   Reference:  N. J. Higham, Computing the nearest correlation
%   matrix---A problem from finance. IMA J. Numer. Anal.,
%   22(3):329-343, 2002.
% 
% Note that the Newton Method algorithm is faster, but this alternating
% projection algorithm is simpler to implement and good enough for this use
% case.

% Sanity check.
assert(isequal(A, A'), 'A must be symmetric.');

% Maximum number of iterations.
max_iter = 512;

% Numerical tolerance for convergence.
tol_con = 16 * eps;

% Numerical tolerance for eigenvalues.
tol_eig = 32 * eps;

% Size of square matrix.
n = length(A);

% Compute initial conditions.
X = A; Y = A;
iter = 0;
rel_diffX = inf; rel_diffY = inf; rel_diffXY = inf;
dS = zeros(size(A));

% Iterate until difference between projections is within tolerance.
while max([rel_diffX rel_diffY rel_diffXY]) > tol_con
    % Keep track of difference between iterations.
    R = Y - dS;
    
    % Projection with positive eigenvalues.
    Xold = X;
    [V,D] = eig(R);
    X = V*diag(max(diag(D), tol_eig))*V';
    X = (X+X')/2; % ensure symmetry
    
    % Keep track of differnce between interations.
    dS = X - R;
    
    % Projection with unit diagonal.
    Yold = Y;
    Y = X;
    Y(1:n+1:n^2) = 1;
    
    % Check for convergence and iterate.
    rel_diffX = norm(X-Xold,'fro')/norm(X,'fro');
    rel_diffY = norm(Y-Yold,'fro')/norm(Y,'fro');
    rel_diffXY = norm(Y-X,'fro')/norm(Y,'fro');
    iter = iter + 1;
    
    % Abort if we have exceeded maximum # of iterations.
    if iter > max_iter
        %warning('nearcorr exceeded maximum number of iterations');
        break;
    end
end