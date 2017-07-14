function A = prepare_graph(A)
% PREPARE_GRAPH Remove diagonals, edge weights, and make sure it's symmetric
%
% A = prepare_graph(A);

A = A|A';
A = A - diag(diag(A));
