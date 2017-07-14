function [result bounds time nmatvec] = commute_pairwise(B, i, j, tol, l, lambda1)
% COMMUTE_PAIRWISE Compute bounds on the commute time between nodes
%
% [result bounds time] = commute_pairwise(A,i,j,tol,maxiter,lambda1)
%
% Example:
%
% See also LAPLACIAN

% History
% :2010-06-23: David's version based on Pooya's old codes

B = laplacian(B);

start_time = tic;
[result bounds nmatvec] = comprob1(B, i, j, tol, l, lambda1);
time = toc(start_time);


function [result bounds nmatvec] = comprob1(A, i, j, tol, l, lambda1)

[ro co] = size(A);
e_i = zeros(ro,1);
e_j = zeros(ro,1);

e_i(i,1) = 1;
e_j(j,1) = 1;

v = e_i - e_j;
u = v / norm(v);

if (nargin < 5)
    l = ceil(sqrt(size(A, 1)));
end

out = zeros(2,l);

[b b1 b2 b3 l] = lapgql(A, u, l, tol, lambda1);

out(1, 1:l) = max([b; b2]) .* 2; %lower bound
out(2, 1:l) = min([b1; b3]) .* 2; % upper bound

bounds = out(:,1:l);
result = (out(1,l) + out(2,l))/2;

nmatvec = l;

function [b b1 b2 b3 j] = lapgql(A, u, l, tol, lambda1, lambda2) %Laplacian Version!!!!

if (nargin < 4)
    tol = 0;
end

if (nargin < 5)
    E = eigs(A, 3, 'sa'); %opts?
    lambda1 = E(2)/1.1;
%    lambda1 = 1e-8; 
end

if (nargin < 6)
    lambda2 = size(A, 1); 
end

if (nargin < 3)
    l = ceil(sqrt(size(A, 1)));
end

row  = size(A, 1);

omega = sparse(l);
omega1 = sparse(l);
omega2 = sparse(l);
omega3 = sparse(l);
gamma = sparse(l);
gamma3 = sparse(l);
c = sparse(l);
d = sparse(l);
d1 = sparse(l);
d2 = sparse(l);
b = sparse(l); 
b1 = sparse(l);
b2 = sparse(l);
b3 = sparse(l);
h = sparse(row,3);
h1 = sparse(row,3);
e = ones(row, 1) / sqrt(row); %%%%

h_1 = 0;
h0 = u;
Au = A * u;
omega = u' * (Au + e * (e' * u));
gamma = norm(Au + e * (e' * u) - omega * u);
b = 1/omega;
d = omega;
c = 1;
d1 = omega - lambda1;
d2 = omega - lambda2;
h(:,2) = (Au + e*(e' * u)  - omega * u )/ gamma;
j = 1;
for j=2:l
  %if (mod(j, 20)==1)
     % j
  %end
  
  Ah = A * h(:,mod(j-1, 3)+1);

  omega(j) = h(:,mod(j-1, 3)+1)' * (Ah + e * (e' * h(:,mod(j-1, 3)+1)));
  if j == 2
    h1(:,3) = (Ah + e * (e' * h(:,2))  - omega(2) * h(:,2)) - gamma(1) * h0;
  else
    h1(:,mod(j, 3)+1) = (Ah + e * (e' * h(:,mod(j-1, 3)+1))  - omega(j) * h(:,mod(j-1, 3)+1)) - gamma(j-1) * h(:,mod(j-2, 3)+1);
  end

  gamma(j) = norm (h1(:,mod(j, 3)+1));

  if gamma(j)==0
      h(:,mod(j, 3)+1) = 0;
  else
      h(:,mod(j, 3)+1) = h1(:,mod(j, 3)+1) / gamma(j);
  end

  % a lower bound b_j of u^TA^-1u by the Gauss quadrature rule
  b(j) = b(j-1) + (gamma(j-1)^2 * c(j-1)^2) / ( (d(j-1) * (omega(j) * d(j-1)-gamma(j-1)^2)) );

  d(j) = omega(j) - gamma(j-1)^2 / ( d(j-1) );

  c(j) = c(j-1) * gamma(j-1) / ( d(j-1));

  d1(j) = omega(j) - lambda1 - gamma(j-1)^2 / (d1(j-1));

  d2(j) = omega(j) - lambda2 - gamma(j-1)^2 / (d2(j-1));

  omega1(j) = lambda1 + gamma(j)^2 / (d1(j));

  omega2(j) = lambda2 + gamma(j)^2 / (d2(j));

  % an upper bound b_j through the Gauss-Radau quadrature rule
  b1(j) = b(j) + (gamma(j)^2 * c(j)^2) / ( (d(j) * (omega1(j) * d(j) - gamma(j)^2)) );

  % a lower bound b_j through the Gauss-Radau quadrature rule
  b2(j) = b(j) + (gamma(j)^2 * c(j)^2) / ( (d(j) * (omega2(j) * d(j) - gamma(j)^2)) );

  omega3(j) = d1(j) * d2(j) / ((d2(j)) - (d1(j))) * (lambda2 / (d1(j)+eps) - lambda1 / (d2(j)));

  gamma3(j) = d1(j) * d2(j) / ((d2(j)) - (d1(j))) * (lambda2 - lambda1);

  % an upper bound b_j through the Gauss-Lobatto rule
  b3(j) = b(j) + (gamma3(j) * c(j)^2) / ((d(j) * (omega3(j) * d(j)-gamma3(j))));  
  
  if (j > 2) && ((relerr(b, j) < tol) && (relerr(b1, j) < tol) && ...
	(relerr(b2, j) < tol) && ...
    (relerr(b3, j) < tol) || (b(j) > b1(j)) || (b(j) > b3(j)) || ...
            (b2(j) > b1(j)) || (b2(j) > b3(j)))
      break;
  end

end

function err = relerr(A, i)
if A(i) ~= 0
    err = abs(A(i)-A(i-1))/abs(A(i));
else
    err = inf;
end