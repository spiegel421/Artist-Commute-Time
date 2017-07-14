function [t,pushstats,x] = F_topk(A,i,varargin)
% TODO write documentation
%
% [t,pushstats,x] = F_topk(A,i,'key',value,...)
%
% 'tol' : stopping tolerance (default 1e-5)
% 'omega' : parameter of richardson (default 1.01)
% 'maxpush' : an integer multiple of the number of nonzeros (default 20)
%
% t = the ith row of the F-measure
%
% pushstats.hist = convergence of pushes
% pushstats.r = residual of push scheme
% pushstats.nverts = number of vertices examined
%
% x should approximate the solution L^+ e_i
%


error(nargchk(2,Inf,nargin,'struct'));

n = size(A,1);
if i<1 || i>n,
    error('F_topk:argumentNotInRange',...
        'the value of i (%i) must be in [1,%i]',i,n);
end

di = full(sum(A(:,i)));

opts = struct('maxpush',20, 'omega', 1.01, 'tol', 1e-6);

optsu = struct(varargin{:});
for fn=fieldnames(opts)'
    f=fn{1};
    if isfield(optsu,f), opts.(f) = optsu.(f); end
end

opts.maxpush = ceil(nnz(A)*opts.maxpush);


%[x,r,hist,nverts] = Lp_push_markov_mex(A,i,omega,tol,maxsteps)

deg = sum(A,2);
[x,r,hist,nverts] = Lp_push_markov_bal_mex(A,i,opts.omega,opts.tol,opts.maxpush);
t = -(deg.*x + (deg(i)).*x);
pushstats.hist = hist;
pushstats.resid = r;
pushstats.nverts = nverts;