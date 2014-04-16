function bestre = csf_er(s, samplex, n, m, k)
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

%% Error reduction
[sosample, in] = sort(samplex);
difre = zeros(1, n);
sampledif = s(in) - re(sosample);
for j = 1 : m
    difre(sosample(j)) = sampledif(j);
    if (j > 1)
        gap = (sampledif(j) - sampledif(j - 1)) / (sosample(j) - sosample(j - 1));
        for q = sosample(j - 1) + 1 : sosample(j)
            difre(q) = difre(q - 1) + gap;
        end
    end
end
if (sosample(1) > 1)
    gap = difre(sosample(1)) - difre(sosample(1) + 1);
    for i = sosample(1) - 1 : -1 : 1
        difre(i) = difre(i + 1) + gap;
    end
end
if (sosample(m) < n)
    gap = difre(sosample(m)) - difre(sosample(m) - 1);
    for i = sosample(m) + 1 : 1 : n
        difre(i) = difre(i - 1) + gap;
    end
end
bestre = re' + difre;