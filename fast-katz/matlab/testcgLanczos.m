function testcgLanczos( n,rtol )

%        testcgLanczos( n,rtol )
% generates a positive-definite matrix A
% a solution x and a rhs vector b = A*x
% and asks cgLanczos to solve A*x = b.

% 22 Oct 2007: First version of testcgLanczos.

rand('state',0)
B = rand(n,n);
A = B'*B;
x = ones(n,1);
b = A*x;

show   = true;
check  = true;
itnlim = 10*n;

matrix = false;

if matrix             % Treat A as an explicit matrix
  [x,istop,itn,Anorm,Acond,rnorm,xnorm,D] = cgLanczos( A,b,show,check,itnlim,rtol );

else                  % Treat A as an operator via Amat.m
  Aop = @(x) Amat(A,x);
  [x,istop,itn,Anorm,Acond,rnorm,xnorm,D] = cgLanczos( Aop,b,show,check,itnlim,rtol );
end

if n<=50
  Ainv  = inv(A);
  Dtrue = diag(Ainv);
  disp(' ')
  disp('     Dtrue      D')
  disp([Dtrue D])
  disp(' ')
  [errmax,imax] = max(abs(Dtrue-D));
  disp('Diagonal with maximum error')
  disp([Dtrue(imax) D(imax)])
  disp('Index:');  disp(imax)
  disp('Error:');  disp(errmax)
  keyboard
end
