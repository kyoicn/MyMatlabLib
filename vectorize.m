%% vectorization

function s = vectorize(c)
s = c(:);
% si = size(c);
% n1 = si(1);
% n2 = si(2);
% s = zeros(1, n1 * n2);
% for i = 1 : n1
%     for j = 1 : n2
%         if (mod(i, 2) == 1)
%             s((i - 1) * n1 + j) = c(i, j);
%         else
%             s(i * n1 - j + 1) = c(i, j);
%         end
%     end
% end