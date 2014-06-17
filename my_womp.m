%OMP Recovery w/ weights
%
% stop: stop condition

function hat_y=my_womp(s, samplex, T, k, w, stop)

si = size(T);                                     %  �۲�����С
m = si(1);
N = si(2);
hat_y = zeros(N, 1);                                 %  ���ع�������(�任��)����                     
r_n = s;
Aug_t = zeros(m, k);
pos_array = zeros(k, 1);

if (nargin < 6)
    stop = 0;
end
if (nargin < 5)
    w = ones(N, 1);
end

wmf = diag(w);
wmp = wmf(samplex, samplex);

for times=1:k
%     product = abs(T'*r_n);
%     product = dot(abs(T'*r_n), w);
    product = abs(T'*wmp*r_n);
    [~, pos]=max(product);                       %  ���ͶӰϵ����Ӧ��λ��
    pos_array(times) = pos;                         %  ��¼���ͶӰϵ����λ��
    Aug_t(:, times) = T(:, pos);
    T(:, pos) = zeros(m, 1);
    part = Aug_t(:, 1 : times);
%     aug_y = inv(part'*part)*part'*s;
%     size(wmp)
%     size(part)
    aug_y = pinv(wmp * part)*(wmp * s);
    r_n = s - part * aug_y;                            %  �в�
    
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
    
hat_y(pos_array)=aug_y;                           %  �ع�������
