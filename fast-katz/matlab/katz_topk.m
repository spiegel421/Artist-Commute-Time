function [x,hist,r,nverts,nedges] = katz_topk(A,alpha,i,varargin)
% KATZ_TOPK Estimate the topk closest nodes in terms of Katz score
%
% k = katz_topk(A,i) returns a good estimate of the $i$th row of the
% Katz matrix, i.e. the Katz score between node i and all other
% nodes in the graph.
%
% k = katz_topk(A,i,'key',value,'key',value,...) changes the default
% settings.  The value keys and their default values are:
%   'tol' - the stopping tolerance [ float { 1e-4 } ]
%   'maxpush' - the maximum number of effective matvecs 
%     [ positive integer { 20 } ].  Formally, the number of pushes
%     performed by the algorithm is 20*n where n is the number of vertices
%     in the graph.
%   'alg' - choice of algorithm to estimate the topk results
%     [ {'balanced'} | 'unbalanced' ]  The balanced algorithm uses
%     the vertex with the largest weighted residual, or balanced residual,
%     where the weighting is the inverse degree.  i.e. we pick vertices
%     with the largest r(i)/d(i), where d(i) is the degree of node i.
%
% [k,hist,x,nverts,nedges] = katz_topk(...) also outputs information on the
% history of the computation after each push:
%   hist(i,1) = the number of edges seen
%   hist(i,2) = the largest residual (or balanced residual if 'alg' = balanced)
%   hist(i,3) = the current sum of the residual elements
%   hist(i,4) = the time taken until the ith step.
%
% This function uses the ability of cgLanczos.m to estimate the diagonal
% elements of the inverse of a matrix.
%
% See also KATZ_PAIRWISE
%
% Example:
%   TODO

% History
% :2011-03-22: Initial coding

error(nargchk(3,Inf,nargin,'struct'));

n = size(A,1);
if i<1 || i>n,
    error('katz_topk:argumentNotInRange',...
        'the value of i (%i) must be in [1,%i]',i,n);
end

if alpha<0
    error('katz_topk:argumentNotInRange',...
        'the value of alpha (%f) must be in [0,Inf]',alpha);
end

opts = struct('tol',1e-4, 'maxpush',20, 'alg','balanced');
optsu = struct(varargin{:});
for fn=fieldnames(opts)'
    f=fn{1};
    if isfield(optsu,f), opts.(f) = optsu.(f); end
end

npush = n*opts.maxpush;

switch opts.alg
    case 'balanced'
        [x,r,hist,nverts] = katz_push_bal_mex(A,alpha,i,opts.tol,npush);
        
    case 'unbalanced'
        [x,r,hist,nverts] = katz_push_mex(A,alpha,i,opts.tol,npush);
        
    otherwise
        error('katz_topk:invalidParameter','alg %s is not a known algorithm',...
            opts.alg);
end

hist = hist';
nedges = hist(end,1);


