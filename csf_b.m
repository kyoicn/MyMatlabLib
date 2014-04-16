function [re, samplex] = csf_b(s, m, k, b)
n = length(s);
origin = s';

randp = randperm(n);
samplex = randp(1:m);

%% function base
a = b(samplex, :);

%% OMP
co = omp(origin(samplex)', a, n, k);

%% WOMP
% d = dis(n);
% co = womp(origin(samplex)', a, n, k, d);

%% BP
% rco = bp(origin(samplex)', a, n);
% co=rco';

re = b*co';