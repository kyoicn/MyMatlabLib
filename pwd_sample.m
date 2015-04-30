% Distributed Pairwise Sampling
% So far only for square cost map
% [samples, realm] = pwd_sample(m, cost, r)
% m: sample size
% cost: cost map
% r: radius map

function [samples, realm, message_count] = pwd_sample(m, cost, r)
dim = size(cost);
n = dim(1) * dim(2);

m0 = m;
if (2*m0 > n) m = n - m; end
p = 2*m / n;
cands = find(binornd(1, p, [1, n]) == 1);

mm = length(cands);
samples = zeros(dim);
paired = samples;

message_count = 0;

cox = mod(cands - 1, dim(1)) + 1;
coy = floor((cands - 1)./dim(1)) + 1;
[tx1 tx2] = meshgrid(cox, cox);
[ty1 ty2] = meshgrid(coy, coy);
dis = (tx2-tx1).^2 + (ty2-ty1).^2 + diag(inf*ones(1, mm));
dis(dis > r) = inf;
for sam = 1 : mm
    if (paired(cox(sam), coy(sam)) == 1) continue; end % previously paired
%     idx1 = hcindex(dim(1), cox(sam), coy(sam));
	peer = find(dis(sam, :)<inf);
	if (~isempty(peer))
        pidx = randi(length(peer));
        message_count = message_count + 2;
% 		idx2 = hcindex(dim(1), cox(peer(1)), coy(peer(1)));
		if (cost(cox(sam), coy(sam)) < cost(cox(peer(pidx)), coy(peer(pidx))))
			if (2*m0 > n) samples(cox(peer(pidx)), coy(peer(pidx))) = 1;
			else samples(cox(sam), coy(sam)) = 1; end
		else
			if (2*m0 > n) samples(cox(sam), coy(sam)) = 1;
			else samples(cox(peer(pidx)), coy(peer(pidx))) = 1; end
        end
        paired(cox(sam), coy(sam)) = 1; paired(cox(peer(pidx)), coy(peer(pidx))) = 1;
		dis([sam peer(pidx)], :) = inf*ones(2, mm);
		dis(:, [sam peer(pidx)]) = inf*ones(mm, 2);
	else
		% no peer
        % flip a coin
        accept = binornd(1, 0.5);
        if (accept == 1)
            samples(cox(sam), coy(sam)) = 1;
            paired(cox(sam), coy(sam)) = 1;
            dis(sam, :) = inf(1, mm);
            dis(:, sam) = inf(mm, 1);
        end
	end
end
if (2*m0 > n)
	samples = ones(dim) - samples;
end
realm = sum(samples(:));
end