function [x,r,hist,nverts] = Lp_push(L,ei,tol,maxsteps)
% Lp_push Implement a push strategy for solving Laplacian systems
%
% Solve Lx=(ei-ej) by exploiting the graph nature of
% L to only access vertex out-degrees.
%
% David F. Gleich
% University of British Columbia, 2010

% History
% :2010-01-29: Initial coding based on commute_push


if ~exist('tol','var') || isempty(tol), tol = 1e-3; end
if ~exist('maxsteps','var') || isempty(maxsteps), maxsteps = 40*size(L,1); end
if ~exist('omega','var') || isempty(omega), omega=1; end

% Algorithm: We store the residual vector in an updatable max-heap
% structure.  Then, we examine the vertex with the largest residual and
% distribute its residual to the neighbors.

% convert graph to csr
[rp,ci,ai] = sparse_to_csr(L);
rho = 1;
ai = ai./rho; % make it diagonally dominant

N=length(rp)-1; % NUMBER of verties
% x is the solution vector
x=zeros(N,1);
% r is the resiudal vector
r=zeros(N,1);
r(ei) = 1/rho;

% circular queue data structure
head = 1; tail = 1;
T=zeros(N,1); L=zeros(N,1);

% setup structure for the residual
nedges = 0;
hist = zeros(maxsteps,4); % nedges, max, sumresid, time
dt = tic; % start the timer
visited = zeros(N,1);
sumresid = 0;
normresid = 0;

% add elements from the current residual to the queue
for w=find(r)'
    sumresid = sumresid + r(w);
    normresid = normresid + r(w)*r(w);
    % add r(w) to the heap
    T(tail) = w;
    L(w) = tail;
    tail = tail+1;
    if tail>N, tail=1; end
end


step = 1;
for step=1:maxsteps
    % get the head of the queue
    v = T(head);
    head = head+1;
    L(v) = 0;
    if head>N, head=1; end
    
    visited(v) = 1;
    
    % Distribut it's rank to everyone else.
    
    % increment our rank
    x(v) = x(v) + omega*r(v);
    % set the value to push out to others and kill our residual
    val = omega*r(v);
    %sumresid = sumresid - omega*r(v); % decrement the residual sum
    
    % for each vertex adjacent to v, push to it!
    for ei=rp(v):rp(v+1)-1            % ei is the edge index
        w=ci(ei);          % w is the target
        
        
        normresid = normresid - r(w)*r(w);
        r(w)=r(w)-val*ai(ei); % increase the residual
        sumresid = sumresid - val*ai(ei); % increase its sum
        normresid = normresid + r(w)*r(w);
        
        
        % check if w is in the queue
        if L(w) == 0,
            % w is not in the queue
            if abs(r(w))*rho > tol,
                % add it
                T(tail) = w;
                L(w) = tail;
                tail = tail+1;
                if tail>N, tail=1; end
            end
        else
            % w is in the queue, ignore it
        end
    end
    
    nedges = nedges + (rp(v+1) - rp(v));
    
    hist(step,1) = nedges;
    if head~=tail, resid = abs(r(T(head)))*rho; else resid = 0; end
    hist(step,2) = resid;
    hist(step,3) = sqrt(normresid);
    hist(step,4) = toc(dt);
    
    if resid < tol
        break;
    end
end

hist = hist(1:step,:); % truncate hist
nverts = sum(visited);