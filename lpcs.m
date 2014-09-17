function x = lpcs(c, A, m, alpha)
%x= lpcs(c, A, m, alpha)
%x: a column vector, the probability of choosing sample i
%c: cost vector
%A: sensing matrix
%m: expected number of samples
%alpha: constraint on column sum
%________________________________
% min  c^T * x
% s.t. sum_x >= m
%      abs(A^T * x) <= alpha
%      x_i <= 1
%________________________________

B = [-ones(1, size(A, 1)); A'; -A'; eye(size(A,1), size(A, 1))];
b = [-m; ones(size(A, 2),1) * alpha; ones(size(A, 2), 1) * alpha; ones(size(A, 1), 1)];
x = linprog(c, B, b);