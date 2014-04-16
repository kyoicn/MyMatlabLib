%% distribution

function r = dis(n)
r = zeros(1, n);

%% random sum = 1
% s = 1;
% for i = 1 : n
%     r(i) = s * rand();
%     s = s - r(i);
% end
% r = sort(r, 'descend');

%% m-uniform
% m = round(n / 3);
% for i = 1 : m
%     r(i) = 1 / m;
% end

%% arithmetic progression
m = 5;
d = (m - 1)/(n - 1) * 2;
for i = 1 : n
    r(n - i + 1) = (i - 1) * d / m / n;
%     r(n - i + 1) = 1 + (i - 1) * d / m / n;
end