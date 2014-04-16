%% function poly re

function re = polyre_1d_dct_ext(s, p)

n = length(s);
origin = s;

[trans, b] = dct4(origin);

j = 0;
for i = 1 : n
    if (binornd(1, p) == 1)
        j = j + 1;
        samplex(j) = i;
    end
end

samplen = length(samplex);

%% CDG
% [cdg_trans, cdg_b] = dct4(origin(samplex));
% m = round(samplen * 0.4);
% 
% for i = 1 : m
%     for j = 1 : samplen
%         aa(i,j) = normrnd(1, 1);
%     end
% end
% 
% cdg_co = omp(aa * origin(samplex)', aa * cdg_b, samplen, m);
% cdg_re = cdg_b * cdg_co';

%% Poly recovery
for i = 1 : samplen
    for j = 1 : n
        a(i,j) = sqrt(2 / n) * cos(pi / n * (j - 0.5) * (samplex(i) - 0.5));
    end
end

%with pre-processing
% q = (orth(a'))';
% t = q * pinv(a);
% s = t * origin(samplex)';
% co = omp(s, q, n, 5);
% co = omp(origin(samplex)', a, n, round(samplen/5));

d = dis(n)
co = womp(origin(samplex)', a, n, round(samplen/4), d);

% rco = bp(origin(samplex)', a, n);
% co=rco';
% rre = zeros(1, n);
% for i = 1 : n
%     for j = 1 : n
%         rre(i) = rre(i) + rco(j) * sqrt(2 / n) * cos(pi / n * (j - 0.5) * (i - 0.5));
%     end
% end

% bd = round(n * 0.1);
% co = [rco(1 : bd); zeros(n - bd, 1)];

% co = gmp(origin(samplex)', a, -100, 100);
re = zeros(1, n);
for i = 1 : n
    for j = 1 : n
        re(i) = re(i) + co(j) * sqrt(2 / n) * cos(pi / n * (j - 0.5) * (i - 0.5));
    end
end

figure;
bar(1:n, trans);
hold on;
plot(1:n, co, 'r*');

figure;
bar(1:n, trans-co');

% figure;
% bar(1:n, rco);

figure;
plot(1:n,origin);
hold on;
plot(1:n,re,'k+');
plot(samplex,origin(samplex),'o');
% plot(1:n,rre,'m*');

figure;
plot(1:n, origin - re);
hold on;
z = zeros(1, n);
plot(samplex, z(samplex), 'o');