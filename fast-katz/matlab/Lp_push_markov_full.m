function [x,r,hist,nverts] = Lp_push(A,ei,tol,maxsteps)
% This function was never finished, it was an idea, but 
% another idea worked better (Lp_push_markov)
return
% Lp_push Implement a push strategy for solving Laplacian systems
%
% The idea here is to solve
% (I - P)x = e - ei
% but we remove row ei, so we are solving
% (I - Pbar)y = e
% Thus, we give everything a residual of 1 initially, but we only
% add neighbors of ei to the queue initially.
%
% David F. Gleich
% University of British Columbia, 2010

% History
% :2010-01-29: Initial coding based on commute_push, added regularization 

% convert graph to csr
[rp,ci,ai] = sparse_to_csr(A);
N=length(rp)-1; % NUMBER of verties

%if ~exist('alpha','var') || isempty(alpha), alpha = 1/sqrt(N); end
if ~exist('tol','var') || isempty(tol), tol = 1e-3; end
if ~exist('maxsteps','var') || isempty(maxsteps), maxsteps = 40*N; end
if ~exist('omega','var') || isempty(omega), omega=1; end

% Algorithm: We store the residual vector in an updatable max-heap
% structure.  Then, we examine the vertex with the largest residual and
% distribute its residual to the neighbors.

source = ei;

%ai = ai./rho; % make it diagonally dominant
d = diff(rp);
id = 1./d;

% x is the solution vector
x=zeros(N,1);
% r is the resiudal vector
r=ones(N,1);
r(ei) = 0;



% heap data structure
n=0; % n size of heap 
T=zeros(N,1); L=zeros(N,1);

% setup structure for the residual
nedges = 0;
hist = zeros(maxsteps,4); % nedges, max, sumresid, time
dt = tic; % start the timer
visited = zeros(N,1);
sumresid = n-1;

% add neighbors of ei to the heap
for ei=rp(source):rp(source+1)-1            % ei is the edge index
    w = ci(ei);
    if w == source, continue; end
    % add r(w) to the heap
    n=n+1; T(n)=w; L(w)=n; k=n;
    % move the element up the heap
    j=k; tj=T(j);
    while j>1,                       % j==1 => element at top of heap
        j2=floor(j/2); tj2=T(j2);    % parent element
        if abs(r(tj2))>abs(r(tj)), break;      % parent is larger, so done
        else                         % parent is smaller, so swap
            T(j2)=tj; L(tj)=j2; T(j)=tj2; L(tj2)=j; j=j2;
        end
    end  
end


step = 1;
for step=1:maxsteps
    v=T(1); ntop=T(n); T(1)=ntop; L(ntop)=1; n=n-1; % pop the head off the heap
    L(v) = 0; % remove v from the heap
    k=1; kt=ntop;                   % move element T(1) down the heap
    while 1,
        i=2*k; 
        if i>n, break; end          % end of heap
        if i==n, it=T(i);           % only one child, so skip
        else                        % pick the largest child
            lc=T(i); rc=T(i+1); it=lc;
            if abs(r(rc))>abs(r(lc)), i=i+1; it=rc; end % right child is larger
        end
        if abs(r(kt))>abs(r(it)), break;     % at correct place, so end
        else T(k)=it; L(it)=k; T(i)=kt; L(kt)=i; k=i; % swap
        end
    end                             % end heap down
    
    % check, this should sort the vector v
    % fprintf('Pushing from %4i val %f\n', v, d(v));
    visited(v) = 1;
    
    % now, v is the element with the largest residual.  Let's 
    % distribute it's rank to everyone else.
    
    % set the value to push out to others and kill our residual
    val = r(v);
    % increment our rank
    x(v) = x(v) + val;
    r(v) = 0;
    sumresid = sumresid - val;

    %a=id(v);
    
    % for each vertex adjacent to v, push to it!
    for ei=rp(v):rp(v+1)-1            % ei is the edge index
        w=ci(ei);          % w is the target
        
        if w == source, continue; end % remove this node


        r(w)=r(w)+val*id(w); % increase the residual
        sumresid = sumresid + val*id(w); % increase its sum
        
        % check if w is in the heap
        k=L(w); onlyup=0; 
        if k==0
            % element not in heap, only move the element up the heap
            n=n+1; T(n)=w; L(w)=n; k=n; kt=w; onlyup=1;
        else kt=T(k);
        end
        % update the heap, move the element down in the heap
        while 1 && ~onlyup,
            i=2*k; 
            if i>n, break; end          % end of heap
            if i==n, it=T(i);           % only one child, so skip
            else                        % pick the largest child
                lc=T(i); rc=T(i+1); it=lc;
                if abs(r(rc))>abs(r(lc)), i=i+1; it=rc; end % right child is larger
            end
            if abs(r(kt))>abs(r(it)), break;      % at correct place, so end
            else T(k)=it; L(it)=k; T(i)=kt; L(kt)=i; k=i; % swap
            end
        end
        % move the element up the heap
        j=k; tj=T(j);
        while j>1,                       % j==1 => element at top of heap
            j2=floor(j/2); tj2=T(j2);    % parent element
            if abs(r(tj2))>abs(r(tj)), break;      % parent is larger, so done
            else                         % parent is smaller, so swap
                T(j2)=tj; L(tj)=j2; T(j)=tj2; L(tj2)=j; j=j2;
            end
        end  
    end
    
    nedges = nedges + (rp(v+1) - rp(v));
    
    hist(step,1) = nedges;
    if n>0, resid = abs(r(T(1))); else resid = 0; end
    hist(step,2) = resid;
    hist(step,3) = sumresid;
    hist(step,4) = toc(dt);
    
    if resid < tol
        break;
    end
end



x = x./diff(rp); % normalize the solution

hist = hist(1:step,:); % truncate hist
nverts = sum(visited);
