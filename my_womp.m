%OMP Recovery
%
% stop: stop condition

function hat_y=my_womp(s, T, N, k, w, stop)

Size = size(T);                                     %  观测矩阵大小
m = Size(1);                                        %  测量
hat_y = zeros(N, 1);                                 %  待重构的谱域(变换域)向量                     
r_n = s;
Aug_t = zeros(m, k);
pos_array = zeros(k, 1);

if (nargin < 6)
    stop = 0;
end
if (nargin < 5)
    w = ones(m, 1);
end

wm = diag(w);

for times=1:k
%     product = abs(T'*r_n);
%     product = dot(abs(T'*r_n), w);
    product = abs(T'*r_n*wm);
    [~, pos]=max(product);                       %  最大投影系数对应的位置
    Aug_t(:, times) = T(:, pos);
    T(:, pos) = zeros(m, 1);
    part = Aug_t(:, 1 : times);
%     aug_y = inv(part'*part)*part'*s;
    aug_y = pinv(wm * part)*(wm * s);
    r_n = s - part * aug_y;                            %  残差
    pos_array(times) = pos;                         %  纪录最大投影系数的位置
    
    if (stop < 0)
        continue;
    elseif (norm(r_n) <= stop)
        break;
    end
end

if (times < k)
    % remove zeros
    pos_array = pos_array(1 : times);
end
    
hat_y(pos_array)=aug_y;                           %  重构的向量
