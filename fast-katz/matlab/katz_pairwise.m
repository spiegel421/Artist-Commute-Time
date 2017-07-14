function [result bounds time nmatvec] = katz_pairwise(B, i, j, lambda1, alpha, tol, l)


if (nargin < 4)
    tol = 0.0001;
end

r = size(B, 1);

C=B*ones(r, 1);

A = speye(r) - alpha * B;

if (nargin < 5)
    l = ceil(sqrt(size(A, 1)));
end

start_time = tic;
[result bounds nmatvec] = katzgql(A, i, j, tol, l, lambda1);
time = toc(start_time);

function [result bounds nmatvec] = katzgql(A, i, j, tol, l, lambda1)

[ro co] = size(A);
e_i = zeros(ro,1);
e_j = zeros(ro,1);

e_i(i,1) = 1;
e_j(j,1) = 1;


v = e_i + e_j;
u = v / norm(v);



%out = zeros(4,l);
out = zeros(2,l); % [lower bound; upper bound]

[b b1 b2 b3 l1] = gql(A, u, l, tol, lambda1);
b_u1 = b;
b1_u1 = b1;
b2_u1 = b2;
b3_u1 = b3;

u = (e_i - e_j) / norm(e_i - e_j);
[b b1 b2 b3 l2] = gql(A, u, l, tol, lambda1);
b_u2 = b;
b1_u2 = b1;
b2_u2 = b2;
b3_u2 = b3;

nmatvec = l1 + l2;

if l1 <= l2
    l = l2;
    nextra = l2 - l1;
    ext=ones(1,nextra);
    % repeat the info until the end
    b_u1 = [b_u1(1:l1) b_u1(l1)*ext];
    b1_u1 = [b1_u1(1:l1) b1_u1(l1)*ext];
    b2_u1 = [b2_u1(1:l1) b2_u1(l1)*ext];
    b3_u1 = [b3_u1(1:l1) b3_u1(l1)*ext];
else
    l = l1;
    nextra = l1 - l2;
    ext=ones(1,nextra);

    b_u2 = [b_u2(1:l2) b_u2(l2)*ext];
    b1_u2 = [b1_u2(1:l2) b1_u2(l2)*ext];
    b2_u2 = [b2_u2(1:l2) b2_u2(l2)*ext];
    b3_u2 = [b3_u2(1:l2) b3_u2(l2)*ext];
end

% Truncate them all to the new length
b_u1 = b_u1(1:l);
b1_u1 = b1_u1(1:l);
b2_u1 = b2_u1(1:l);
b3_u1 = b3_u1(1:l);

b_u2 = b_u2(1:l);
b1_u2 = b1_u2(1:l);
b2_u2 = b2_u2(1:l);
b3_u2 = b3_u2(1:l);

% the max/min here are correct (2011-03-18)
%   notes:
%   we want l <= a-b <= u
%   and we have la <= a <= ua; lb <= b <= ub
%   so we have la - ub <= a - b <= ua - lb
%   or max(la1,la2) - min(ub1,ub2) <= a - b <= min(ua1,ua2) - max(lb1,lb2)
out(1,1:l) = (max([b_u1; b2_u1]) - min([b1_u2; b3_u2])) ./ 2; %lower bound
out(2,1:l) = (min([b1_u1; b3_u1]) - max([b_u2; b2_u2])) ./ 2; % upper bound

bounds = out(:,1:l);
result = (out(1,l) + out(2,l))/2;

function [b b1 b2 b3 j] = gql(A, u, l, tol, lambda1, lambda2)

if (nargin < 4)
    tol = 0;
end

if (nargin < 5)
    E = eigs(A, 3, 'sa'); %opts?
    lambda1 = E(1)/1.1
end

if (nargin < 6)
    lambda2 = norm(A, inf); 
end

if (nargin < 3)
    l = ceil(sqrt(size(A, 1)));
end

row  = size(A, 1);

omega = zeros(1,l);
omega1 = zeros(1,l);
omega2 = zeros(1,l);
omega3 = zeros(1,l);
gamma = zeros(1,l);
gamma3 = zeros(1,l);
c = zeros(1,l);
d = zeros(1,l);
d1 = zeros(1,l);
d2 = zeros(1,l);
b = zeros(1,l); 
b1 = zeros(1,l);
b2 = zeros(1,l);
b3 = zeros(1,l);
h = zeros(row,3);
h1 = zeros(row,3);

h_1 = 0;
h0 = u;

Au = A*u;
omega = u' * Au;
gamma = norm(Au - omega*u);
b = inv(omega);
d = omega;
c = 1;
d1 = omega - lambda1;
d2 = omega - lambda2;
h(:,2) = (Au - omega*u)/gamma;
j = 1;

for j=2:l
 % if (mod(j, 20)==1)
 %     j
 % end
    
  Ah = A * h(:, mod(j-1, 3)+1);
  omega(j) = h(:, mod(j-1, 3)+1)' * Ah;
 
  if j == 2
    h1(:, mod(j, 3)+1) = Ah - omega(j)*h(:, mod(j-1, 3)+1) - gamma(j-1) * h0;
  else
    h1(:, mod(j, 3)+1) = Ah - omega(j)*h(:, mod(j-1, 3)+1) - gamma(j-1) * h(:, mod(j-2, 3)+1);
  end

  gamma(j) = norm (h1(:, mod(j, 3)+1));

  if gamma(j)==0
      h(:, mod(j, 3)+1) = 0;
  else
      h(:, mod(j, 3)+1) = h1(:, mod(j, 3)+1) / gamma(j);
  end

  % a lower bound b_j of u^TA^-1u by the Gauss quadrature rule
  b(j) = b(j-1) + (gamma(j-1)^2 * c(j-1)^2) / ( (d(j-1) * (omega(j) * d(j-1)-gamma(j-1)^2))+eps );

  d(j) = omega(j) - gamma(j-1)^2 / ( d(j-1)+eps );

  c(j) = c(j-1) * gamma(j-1) / ( d(j-1)+eps);

  d1(j) = omega(j) - lambda1 - gamma(j-1)^2 / (d1(j-1)+eps);

  d2(j) = omega(j) - lambda2 - gamma(j-1)^2 / (d2(j-1)+eps);

  omega1(j) = lambda1 + gamma(j)^2 / (d1(j)+eps);

  omega2(j) = lambda2 + gamma(j)^2 / (d2(j)+eps);

  % an upper bound b_j through the Gauss-Radau quadrature rule
  b1(j) = b(j) + (gamma(j)^2 * c(j)^2) / ( (d(j) * (omega1(j) * d(j) - gamma(j)^2)) + eps);

  % a lower bound b_j through the Gauss-Radau quadrature rule
  b2(j) = b(j) + (gamma(j)^2 * c(j)^2) / ( (d(j) * (omega2(j) * d(j) - gamma(j)^2)) + eps);

  omega3(j) = d1(j) * d2(j) / ((d2(j)+eps) - (d1(j))) * (lambda2 / (d1(j)+eps) - lambda1 / (d2(j)+eps));

  gamma3(j) = d1(j) * d2(j) / ((d2(j)+eps) - (d1(j))) * (lambda2 - lambda1);

  % an upper bound b_j through the Gauss-Lobatto rule
  b3(j) = b(j) + (gamma3(j) * c(j)^2) / ((d(j) * (omega3(j) * d(j)-gamma3(j)))+eps);  
  
  %if (j>2) & (relerr(b(j-1:j)-b1(j-1:j)+b2(j-1:j)-b3(j-1:j), 2) < tol*2)
  if (j > 2) && ...
    (relerr(b, j) < tol) && (relerr(b1, j) < tol) && ...
    (relerr(b2, j) < tol) && (relerr(b3, j) < tol)
      break;
  end
end



function err = relerr(A, i)
if A(i) ~= 0
    err = abs(A(i)-A(i-1))/abs(A(i));
else
    err = inf;
end