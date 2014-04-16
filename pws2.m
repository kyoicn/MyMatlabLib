%% pairwise sampling

function samplex = pws2(n, m)
if (m > n) 
    disp('m > n');
    samplex = [];
    return;
end
m0 = m;
if (2*m0 > n) m = n - m; end
% tic;
samplex = zeros(1, m);
cands = randperm(n, 2*m);
dis = zeros(1, m*(2*m-1));
count = 0;
for i = 2 : 2*m
    for j = 1 : i - 1
        count = count + 1;
        dis(count) = abs(cands(i) - cands(j));
    end
end
% dis = abs(repmat(cands, 2*m, 1) - repmat(cands', 1, 2*m)) + diag(inf * ones(1, 2*m));
[~, sin] = sort(dis);
taken = zeros(1, 2*m);
% t1 = toc;
% cox = mod(cands, 32);
% coy = ceil(cands./32);
% [tx1 tx2] = meshgrid(cox, cox);
% [ty1 ty2] = meshgrid(coy, coy);
% dis = (tx2-tx1).^2 + (ty2-ty1).^2 + diag(inf*ones(1, 2*m));
% tic;
sam = 0;
for i = 1 : m*(2*m-1)
    in = sin(i);
    group = floor((sqrt(1+8*in) - 1) / 2);
    offset = in - group*(group + 1) / 2;
    if (offset == 0)
        idx1 = group;
        idx2 = group + 1;
    else
        idx1 = offset;
        idx2 = group + 2;
    end
%     fprintf('i=%d, in=%d, g=%d, of=%d, i1=%d, i2=%d:', i, in, group, offset, idx1, idx2);
    if (taken(idx1) || taken(idx2)) continue;
    else
%         fprintf('yes\n');
        sam = sam + 1;
        if (binornd(1, 0.5) == 1) samplex(sam) = cands(idx1);
        else samplex(sam) = cands(idx2); end
        if (sam == m) break; end
    end
    taken([idx1 idx2]) = 1;
end
% t2 = toc;

if (2*m0 > n)
    bif = ones(1,n);
    bif(samplex) = 0;
    samplex = find(bif==1);
end
% m = m0;
% fprintf('%f, %f\n', t1, t2);
end