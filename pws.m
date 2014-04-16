%% pairwise sampling

function samplex = pws(n, m)
if (m > n) 
    disp('m > n');
    samplex = [];
    return;
end
m0 = m;
if (2*m0 > n) m = n - m; end
samplex = zeros(1, m);
cands = sort(randperm(n, 2*m));

% decision = binornd(1, 0.5, 1, m);
decision = zeros(1, m);
decision(2:2:m) = 1;
for sam = 1 : m
    samplex(sam) = cands(sam*2 - decision(sam));
end

if (2*m0 > n)
    bif = ones(1,n);
    bif(samplex) = 0;
    samplex = find(bif==1);
end
end