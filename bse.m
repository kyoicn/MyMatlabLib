%% binary sparse enumeration

function re = bse(obs, matrix, n, k)
re = zeros(n, 1);
for i = 1:k
    