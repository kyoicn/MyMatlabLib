function [bestre, samplex] = csf_er_exp(s, m, k)
n = length(s);
origin = s';

% [trans, b] = dct4(origin);

randp = randperm(n);
samplex = randp(1:m);

[st, dctb] = dct4(s');

%% Func recovery
a = dctb(samplex, :);

%% OMP with pre-processing
% q = (orth(a'))';
% t = q * pinv(a);
% s = t * origin(samplex)';
% co = omp(s, q, n, k);

%% OMP
co = omp(origin(samplex)', a, n, k);

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

%% Error reduction
[sosample, insample] = sort(samplex,2);
difre = zeros(1, n);
sampledif = s(sosample) - re(sosample);
for j = 1 : m
    difre(sosample(j)) = sampledif(j);
    if (j > 1)
        gap = (sampledif(j) - sampledif(j - 1)) / (sosample(j) - sosample(j - 1));
        for q = sosample(j - 1) + 1 : sosample(j)
            difre(q) = difre(q - 1) + gap;
        end
    end
end

bestre = re' + difre;