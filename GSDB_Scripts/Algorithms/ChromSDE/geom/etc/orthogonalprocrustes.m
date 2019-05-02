% orthogonalprocurstes Solve the orthogonal Procrustes problem.
%
% Solves the orthogonal Procrustes problem:
%            Find U such that U^TU = I, and U minimizes ||A-UB||_F.
%
% U = orthogonalprocrustes(A,B)

% Created       :  6 Jul 07
% Last modified :  6 Jul 07

function U = orthogonalprocrustes(A,B)

[W,S,Z] = svd(A*B');
U = W * Z';
