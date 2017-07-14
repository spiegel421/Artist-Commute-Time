function [r,d,hist,nverts] = katz_push(A,a,b,tol,maxsteps,omega)
% KATZ_PUSH Implement a push strategy for Katz's ranking
%
% [x,r,hist,nverts] = katz_push(A,alpha,b,tol,maxsteps,omega) specifies
%   A - the adjacency matrix for the graph
%   alpha - the value of alpha in the Katz matrix
%   b - the non-zeros in the right hand side
%   tol - a stopping tolerance (default = 1e-7);
%   maxsteps - a limit on the total number of steps (default = 10*size(A,1))
%   omega - the acceleration parameter in richardson, only omega ~ 1 works
% and returns
%   x - the Katz ranking vector given by the solution of
%     (I-aA)x = b
%   r - the residual vector 
%   hist - a history of the residuals at each step
%   nverts - the total number of unique vertices explored
%   
% Solve (I-aA)x=b by exploiting the graph nature of
% A to only access vertex out-degrees.
%
% If b has length==N, then it is a full right hand side, otherwise, it 
% is a list of vertices and we set b = 1 for all these vertices and
% 0 elsewhere.

% David F. Gleich
% University of British Columbia, 2010

% History
% -------
% :2010-01-14: Initial coding with elements of the heap from Dijkstra's
% :2010-01-29: Added documentation


if ~exist('tol','var') || isempty(tol), tol = 1e-7; end
if ~exist('maxsteps','var') || isempty(maxsteps), maxsteps = 10*size(A,1); end
if ~exist('omega','var') || isempty(omega), omega=1; end

% Algorithm: We store the residual vector in an updatable max-heap
% structure.  Then, we examine the vertex with the largest residual and
% distribute its residual to the neighbors.

% convert graph to csr
[rp,ci] = sparse_to_csr(A);

N=length(rp)-1; % NUMBER of verties
% r is the ranking vector
r=zeros(N,1);
% d is the residual vector
if length(b) == N
    d=b(:); 
else
    % otherwise, this is a list of vertices
    d = zeros(N,1);
    d(b) = 1;
end

% heap data structure
n=0; % n size of heap 
T=zeros(N,1); L=zeros(N,1);

% setup structure for the residual
nedges = 0;
hist = zeros(maxsteps,4); % nedges, max, sumresid, time
dt = tic; % start the timer
visited = zeros(N,1);
sumresid = 0;

% add elements from the current residual to the heap
for w=find(d)'
    sumresid = sumresid + d(w);
    % add d(w) to the heap
    n=n+1; T(n)=w; L(w)=n; k=n;
    % move the element up the heap
    j=k; tj=T(j);
    while j>1,                       % j==1 => element at top of heap
        j2=floor(j/2); tj2=T(j2);    % parent element
        if d(tj2)>d(tj), break;      % parent is larger, so done
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
            if abs(d(rc))>abs(d(lc)), i=i+1; it=rc; end % right child is larger
        end
        if abs(d(kt))>abs(d(it)), break;     % at correct place, so end
        else T(k)=it; L(it)=k; T(i)=kt; L(kt)=i; k=i; % swap
        end
    end                             % end heap down
    
    % check, this should sort the vector v
    % fprintf('Pushing from %4i val %f\n', v, d(v));
    visited(v) = 1;
    
    % now, v is the element with the largest residual.  Let's 
    % distribut it's rank to everyone else.
    
    % increment our rank
    r(v) = r(v) + omega*d(v);
    % set the value to push out to others and kill our residual
    val = a*omega*d(v);
    sumresid = sumresid - omega*d(v); % decrement the residual sum
    d(v) = d(v) - omega*d(v);
    
    
    
    % for each vertex adjacent to v, push to it!
    for ei=rp(v):rp(v+1)-1            % ei is the edge index
        w=ci(ei);          % w is the target
        
        d(w)=d(w)+val; % increase the residual
        sumresid = sumresid + val; % increase its sum
        
        %fprintf('%i -> %i with value %f  (r[w]=%f, sumresid=%f, n=%i)\n', ...
        %            v, w, val, d(w), sumresid, n);
        
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
                if abs(d(rc))>abs(d(lc)), i=i+1; it=rc; end % right child is larger
            end
            if abs(d(kt))>abs(d(it)), break;      % at correct place, so end
            else T(k)=it; L(it)=k; T(i)=kt; L(kt)=i; k=i; % swap
            end
        end
        % move the element up the heap
        j=k; tj=T(j);
        while j>1,                       % j==1 => element at top of heap
            j2=floor(j/2); tj2=T(j2);    % parent element
            if abs(d(tj2))>abs(d(tj)), break;      % parent is larger, so done
            else                         % parent is smaller, so swap
                T(j2)=tj; L(tj)=j2; T(j)=tj2; L(tj2)=j; j=j2;
            end
        end  
    end
    
    nedges = nedges + (rp(v+1) - rp(v));
    
    hist(step,1) = nedges;
    if n>0, resid = abs(d(T(1))); else resid = 0; end
    hist(step,2) = resid;
    hist(step,3) = sumresid;
    hist(step,4) = toc(dt);
    
    if sumresid>n*1./tol
        % failure!
        warning('katz_push:residualTooLarge',...
            'after %i steps, the residual sum is %g -- halting iteration',...
            step, sumresid);
        break
    end
    
    if resid < tol
        break;
    end
end

hist = hist(1:step,:); % truncate hist
nverts = sum(visited);