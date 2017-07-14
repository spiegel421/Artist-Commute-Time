function y = Amat(A,x)

%        y = Amat(A,x);
% assumes A is a matrix and allows it to be treated
% as an operator for forming A*x.
%   Aop = @(x) Amat(A,x);
% allows Aop to be passed as a function handle
% to a routine that is expecting an operator.

% 22 Aug 2007: First version of Amat to go with cgLanzcos.

  y = A*x;
  return
