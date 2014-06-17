function [bestre co] = my_csf(s, samplex, base, k, method, er, w)
% [bestre co] = my_csf(s, samplex, base, k, method, er)

a = base(samplex, :);
si = size(base);

n = si(1);
m = length(samplex);
q = si(2);
co = zeros(n, 1);

if (nargin < 7)
    w = ones(size(n, 1));
end

%% Methods
if (strcmp(method, 'OMP') == 1)
    co = my_omp(s, a, q, k, 0);
elseif (strcmp(method, 'WOMP') == 1)
    co = my_womp(s, a, q, k, w, 0);
elseif (strcmp(method, 'BP') == 1)
    co = bp(s, a, q);
end

% OMP with pre-processing
% q = (orth(a'))';
% t = q * pinv(a);
% s = t * origin(samplex)';
% co = omp(s, q, n, k);

% WOMP
% d = dis(n);
% co = womp(origin(samplex)', a, n, k, d);

re = base * co;
bestre = re;

%% Error reduction
if (er == 1)
    [sosample, in] = sort(samplex);
%     sosample = unique(sosample);
%     in = unique(in);
%     difre = zeros(1, n);
    sampledif = s(in) - re(sosample);
    difre = interp1(sosample, sampledif, 1:n, 'linear', 'extrap');
%     for j = 1 : m
%         difre(sosample(j)) = sampledif(j);
%         if (j > 1)
%             gap = (sampledif(j) - sampledif(j - 1)) / (sosample(j) - sosample(j - 1));
%             for q = sosample(j - 1) + 1 : sosample(j)
%                 difre(q) = difre(q - 1) + gap;
%             end
%         end
%     end
%     if (sosample(1) > 1)
%         gap = difre(sosample(1)) - difre(sosample(1) + 1);
%         for i = sosample(1) - 1 : -1 : 1
%             difre(i) = difre(i + 1) + gap;
%         end
%     end
%     if (sosample(m) < n)
%         gap = difre(sosample(m)) - difre(sosample(m) - 1);
%         for i = sosample(m) + 1 : 1 : n
%             difre(i) = difre(i - 1) + gap;
%         end
%     end
    bestre = re + difre';
end