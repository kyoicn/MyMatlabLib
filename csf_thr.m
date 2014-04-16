function re = csf_thr(s, m, k, thr)
n = length(s);
origin = s';
% [trans, b] = dct4(origin);
part = 0.5;
m1 = floor(m * part);
m2 = m - m1;

%% phase-1
randp = randperm(n);
samplex = randp(1:m1);

samplen = m1;

%% Func recovery
a=zeros(m1,n);
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

%% Phase-1 Recovery
r = zeros(1, n);
for i = 1 : n
    for j = 1 : n
        r(i) = r(i) + co(j) * sqrt(2 / n) * cos(pi / n * (j - 0.5) * (i - 0.5));
    end
end

%% phase-2

%% Concentration above thr
% exn = abovethr(r, thr); % number of points expected to be above thr
% part = 0.5;
% r1 = floor(m * part);
% r2 = m - r1;
% samplex = zeros(1, m);
% 
% if (exn > r1 && exn < n - m + r1)
%     randp = randperm(n);
%     i = 1; j = 1;
%     while (i < r1 + 1)
%         if (r(randp(j)) > thr)
%             samplex(i) = randp(j);
%             i = i +1;
%         end
%         j = j + 1;
%     end
%     randp = randperm(n);
%     i = 1; j = 1;
%     while (i < r2 + 1)
%         if (r(randp(j)) <= thr)
%             samplex(r1 + i) = randp(j);
%             i = i +1;
%         end
%         j = j + 1;
%     end
% elseif (exn <= r1)
%     randp = randperm(n);
%     i = 1; j = 1; c = 0;
%     while (i < m + 1)
%         if (r(randp(j)) > thr)
%             samplex(i) = randp(j);
%             i = i +1;
%         else
%             if (c < m - r1)
%                 samplex(i) = randp(j);
%                 i = i + 1;
%             end
%         end
%         j = j + 1;
%     end
% else
%     randp = randperm(n);
%     i = 1; j = 1; c = 0;
%     while (i < m + 1)
%         if (r(randp(j)) <= thr)
%             samplex(i) = randp(j);
%             i = i +1;
%         else
%             if (c < m + exn - n)
%                 samplex(i) = randp(j);
%                 i = i + 1;
%             end
%         end
%         j = j + 1;
%     end
% end

%% Boundary detection
tmpsample = samplex;
samplex = keypoints(r, m, 'GENERAL');
% p = m2 + 1;
% for i = 1 : m1
%     flag = 0;
%     for j = 1 : m2
%         if (tmpsample(i) == samplex(j))
%             flag = 1;
%             break;
%         end
%     end
%     if (flag == 0)
%         samplex(p) = tmpsample(i);
%         p = p + 1;
%     end
% end
realm = length(samplex);

%% Phase-2 matrix generation
a=zeros(realm,n);
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

%% Phase-2 Recovery
re = zeros(1, n);
for i = 1 : n
    for j = 1 : n
        re(i) = re(i) + co(j) * sqrt(2 / n) * cos(pi / n * (j - 0.5) * (i - 0.5));
    end
end
re = re';