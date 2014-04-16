function re = csf(s, samplex, n, m, k)
[st, dctb] = dct4(zeros(n, 1));
a = dctb(samplex, :);

%% OMP with pre-processing
% q = (orth(a'))';
% t = q * pinv(a);
% s = t * origin(samplex)';
% co = omp(s, q, n, k);

%% OMP
co = omp(s, a, n, k);

%% WOMP
% d = dis(n);
% co = womp(origin(samplex)', a, n, k, d);

%% BP
% rco = bp(origin(samplex)', a, n);
% co=rco';
% 
% rre = zeros(1, n);
% for i = 1 : n
%     for j = 1 : n
%         rre(i) = rre(i) + rco(j) * sqrt(2 / n) * cos(pi / n * (j - 0.5) * (i - 0.5));
%     end
% end
% 
% bd = round(n * 0.1);
% co = [rco(1 : bd); zeros(n - bd, 1)];
% 
% co = gmp(origin(samplex)', a, -100, 100);
re = dctb * co';