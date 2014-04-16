%% function poly re

function re = polyre_1d_dct(s,m,k)

n = length(s);
origin = s';

% [trans, b] = dct4(origin);

randp = randperm(n);
samplex = randp(1:m);

samplen = m;

%% Func recovery
a=zeros(m,n);
for i = 1 : samplen
    for j = 1 : n
        a(i,j) = sqrt(2 / n) * cos(pi / n * (j - 0.5) * (samplex(i) - 0.5));
    end
end

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
re=re';
% figure;
% bar(1:n, trans);
% hold on;
% plot(1:n, co, 'r*');
% 
% figure;
% bar(1:n, trans-co');
% 
% % figure;
% % bar(1:n, rco);
% 
% figure;
% plot(1:n,origin);
% hold on;
% plot(1:n,re,'k+');
% plot(samplex,origin(samplex),'o');
% % plot(1:n,rre,'m*');
% 
% figure;
% plot(1:n, origin - re);
% hold on;
% z = zeros(1, n);
% plot(samplex, z(samplex), 'o');