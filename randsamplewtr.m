%% weighted sampling w/out replacement

function idx = randsamplewtr(n, m, w)
samplex = zeros(1, n);
chosen = 0;
mm = m;
while (chosen < m)
    tx = unique(randsample(n, mm, true, w));
    tn = length(tx);
    samplex(tx) = ones(1, tn);
    chosen = chosen + tn;
    mm = m - chosen;
    w(tx) = zeros(1, tn);
end
idx = find(samplex==1);

% if (size(w,1) > size(w,2)) w = w'; end
% u = rand(1, n) .^ (1 ./ w);
% [~, in] = sort(u, 'descend');
% idx = in(1:m);