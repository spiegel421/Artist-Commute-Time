function [c,d,reshist,x] = commute_topk(A,i,varargin)
% COMMUTE_TOPK Estimate the topk closest nodes in terms of commute time
%
% c = commute_topk(A,i) returns a course estimate of the $i$th row of the
% commute time matrix, i.e. the commute times between node i and all other
% nodes in the graph.
%
% c = commute_topk(A,i,'key',value,'key',value,...) changes the default
% settings.  The value keys and their default values are:
%   'tol' - the stopping tolerance [ float { 1e-4 } ]
%   'maxit' - the maximum number of iterations [ positive integer { 1000 } ]
%   'diag' - a prescribed vector for the diagonal elements of the inverse
%     [ a length n vector { computed by cgLanczos } ]
%   'reorth' - the number of vectors to use for reorthogonalization in
%     the cgLanczos routine [ a non-negative integer { 10 } ]
%
% [c,d,reshist,x] = commute_topk(...) also returns the estimate of 
% diagonal, the residual history, along with the solution of the 
% pseudo-inverse system L^+ x = e_i, which is used to estimate 
% the top-k commute times.
%
% This function uses the ability of cgLanczos.m to estimate the diagonal
% elements of the inverse of a matrix.
%
% See also LAPLACIAN, CGLANCZOS

% History
% :2011-03-10: Added key/value option

error(nargchk(2,Inf,nargin,'struct'));

n = size(A,1);
if i<1 || i>n,
    error('commute_topk:argumentNotInRange',...
        'the value of i (%i) must be in [1,%i]',i,n);
end

opts = struct('tol',1e-6, 'maxit',1000, 'diag',[], 'reorth', 10);
optsu = struct(varargin{:});
for fn=fieldnames(opts)'
    f=fn{1};
    if isfield(optsu,f), opts.(f) = optsu.(f); end
end

% construct the preconditioner
deg = sum(A,2);
d2 = sqrt(deg);
id2 = 1./d2;

L = laplacian(A);

Ln = @(x) (id2.*(L*(id2.*x) + mean(id2.*x)));
b = zeros(n,1);
b(i) = 1;



[x,flag,iter,Anorm,Acond,rnorm,xnorm,reshist,d] = ...
    cgLanczos(Ln,id2.*b,0,false,opts.maxit,opts.tol,opts.reorth);

if flag==6
    error('commute_topk:cgLanczosError',...
        'L not positive definite, (is A symmetric?)');
elseif flag == 4
    warning('commute_topk:didNotConverge',...
        'The iteration did not converge to rtol=%g in %i iterations',...
        tol, maxit);
elseif flag ==3
    warning('commute_topk:cgLanczosError',...
        'L is highly ill-conditioned, iteration stopped early.');
elseif flag == 2
    warning('commute_topk:didNotFullyConverge',...
        'The iteration converged to machine precision %g, but not rtol=%g',...
        rnorm, tol);
end
   

if ~isempty(opts.diag)
    d = opts.diag;
else
    % adjust d
    d = id2.*d.*id2;
    d = d - 1/n;
end

% adjust x
x = id2.*x - 1/n;

c = d+x(i)-2*x;
c(i) = 0.; % fix this one at zero
