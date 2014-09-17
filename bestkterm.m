% return best k-term approximation of a given vector

function [vk, vknorm, err] = bestkterm(v, k)
n = length(v);
vk = zeros(size(v));
[~, pos] = sort(abs(v), 'descend');
vk(pos(1:k)) = v(pos(1:k));
vknorm = norm(vk);
err.err2 = norm(v-vk);
err.err1 = norm(v-vk, 1);
end