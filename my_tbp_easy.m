function [y_recover x_recover] = my_tbp_easy(y, sigma, thr, A, eps)

A_full = A;

mask = find(~isnan(y));
y = y(mask);
sigma = sigma(mask);
A = A(mask, :);

keep = find(sigma(:) <= quantile(sigma, thr));

x0 = A(keep,:)'*inv(A(keep,:)*A(keep,:)')*y(keep);
epsilon = sqrt(sigma(keep)' * sigma(keep)) * eps;
x_TBP = l1qc_logbarrier(x0, A(keep, :), [], y(keep), epsilon, 1e-3);
y_TBP = A_full * x_TBP;
fprintf(1,'TBP number of nonzero weights: %d\n',sum(x_TBP~=0));
x_recover = x_TBP;
y_recover = y_TBP;

