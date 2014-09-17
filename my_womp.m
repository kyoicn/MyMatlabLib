%OMP Recovery w/ weights
%
% stop: stop condition

function [hat_y, pos_array, er] = my_womp(s, T, k, w, stop)

[m N] = size(T);
hat_y = zeros(N, 1);                                 %  待重构的谱域(变换域)向量                     
r_n = s;
Aug_t = zeros(m, k);
pos_array = zeros(k, 1);
er = zeros(k, 1);

if (nargin < 4) w = ones(m, 1); end
if (nargin < 5) stop = 0; end

wm = diag(w);

for times=1:k
%     product = abs(T'*r_n);
%     product = dot(abs(T'*r_n), w);
    product = abs(T'*wm*r_n);
    [~, pos]=max(product);                       %  最大投影系数对应的位置
    pos_array(times) = pos;                         %  纪录最大投影系数的位置
    Aug_t(:, times) = T(:, pos);
    T(:, pos) = zeros(m, 1);
    part = Aug_t(:, 1 : times);
    aug_y = inv(part'*part)*part'*s;
%     aug_y = pinv(wm * part)*(wm * s);
    r_n = s - part * aug_y;
    er(times) = norm(r_n);
    
    if (stop < 0)
        continue;
    elseif (norm(r_n) <= stop)
        break;
    end
end

if (times < k)
    % remove zeros
    pos_array = pos_array(1 : times);
    er = er(1:times);
end
    
hat_y(pos_array)=aug_y;
