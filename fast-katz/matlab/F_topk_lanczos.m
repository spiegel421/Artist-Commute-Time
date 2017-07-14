function [f,d,reshist,x] = F_topk_lanczos(A,i,varargin)
% TODO write documentation
%
% [f,d,reshist,x] = F_topk_lanczos(A,i,'key',value,...)
%
% 'tol' : stopping tolerance (default 1e-5)
%
%  = the ith row of the F-measure
% 
%
% pushstats.hist = convergence of pushes
% pushstats.r = residual of push scheme
% pushstats.nverts = number of vertices examined
%

% error(nargchk(2,Inf,nargin,'struct'));
% 
% n = size(A,1);
% if i<1 || i>n,
%     error('F_topk:argumentNotInRange',...
%         'the value of i (%i) must be in [1,%i]',i,n);
% end
% 
% di = full(sum(A(:,i)));
% 
% opts = struct('maxpush',20, 'omega', 1.01, 'tol', 1e-6);
% 
% optsu = struct(varargin{:});
% for fn=fieldnames(opts)'
%     f=fn{1};
%     if isfield(optsu,f), opts.(f) = optsu.(f); end
% end
% 
% opts.maxpush = ceil(nnz(A)*opts.maxpush);
% 
% 
% %[x,r,hist,nverts] = Lp_push_markov_mex(A,i,omega,tol,maxsteps)
% 
% deg = sum(A,2);
% [x,r,hist,nverts] = Lp_push_markov_bal_mex(A,i,opts.omega,opts.tol,opts.maxpush);
% t = -(deg.*x + (deg(i)).*x);
% pushstats.hist = hist;
% pushstats.resid = r;
% pushstats.nverts = nverts;


error(nargchk(2,Inf,nargin,'struct'));

n = size(A,1);
if i<1 || i>n,
    error('commute_topk:argumentNotInRange',...
        'the value of i (%i) must be in [1,%i]',i,n);
end

opts = struct('tol',1e-6, 'maxit',1000, 'reorth', 10);
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
    error('F_topk_lanczos:cgLanczosError',...
        'L not positive definite, (is A symmetric?)');
elseif flag == 4
    warning('F_topk_lanczos:didNotConverge',...
        'The iteration did not converge to rtol=%g in %i iterations',...
        tol, maxit);
elseif flag ==3
    warning('F_topk_lanczos:cgLanczosError',...
        'L is highly ill-conditioned, iteration stopped early.');
elseif flag == 2
    warning('F_topk_lanczos:didNotFullyConverge',...
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
x = id2.*x; % don't subtract off the mean!
assert(min(x) >= 0);

f = -(deg.*x + (deg(i)).*x);

c = d+x(i)-2*x;
c(i) = 0.; % fix this one at zero
