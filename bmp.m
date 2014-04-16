%% binary matching pursuit

function r=bmp(obs, matrix, n, k)
si = size(obs);
m = si(1);
res = obs;
r = [];
layer = 1;
c = 1;
opt = zeros(k, n);
tm = matrix;
pt = zeros(1, k);
sel = zeros(1, k);

while (1)
    if (pt(layer) == 0) %just enter a new layer
        cor = zeros(1, n);
        for i = 1 : n
            cor(i) = res' * tm(:, i);
        end
        mv = max(cor);
        if (mv == 0) %over
            r(c, :) = sel;
            c = c + 1;
            layer = layer - 1;
            if (layer == 0)
                break;
            end
            tm(:, opt(layer, pt(layer) - 1)) = matrix(:, opt(layer, pt(layer) - 1));
            res = mod(res + tm(:, opt(layer, pt(layer) - 1)), 2);
        else
            tp = 1;
            for i = 1 : n
                if (cor(i) == mv)
                    opt(layer, tp) = i;
                    tp = tp + 1;
                end
            end
            pt(layer) = 1;
        end
    elseif (opt(layer, pt(layer)) ~= 0)
        sel(layer) = opt(layer, pt(layer));
        if (layer == k) %over
            r(c, :) = sel;
            c = c + 1;
            pt(layer) = pt(layer) + 1;
        else
            res = mod(res + tm(:, sel(layer)), 2);
            tm(:, sel(layer)) = zeros(m, 1);
            pt(layer) = pt(layer) + 1;
            layer = layer + 1;
        end
    elseif (opt(layer, pt(layer)) == 0) %layer over
        sel(layer) = 0;
        opt(layer, :) = zeros(1, n);
        pt(layer) = 0;
        layer = layer - 1;
        if (layer == 0)
            break;
        end
        tm(:, opt(layer, pt(layer) - 1)) = matrix(:, opt(layer, pt(layer) - 1));
        res = mod(res + tm(:, opt(layer, pt(layer) - 1)), 2);
    end
end
% r = zeros(n, 1);
% res = obs;
% thr = 5;
% si = size(obs);
% m = si(1);
% c = 0;
% pos = [];
% while(c < k && hweight(res) >= thr)
%     c = c + 1;
%     for i = 1 : n
%         if (any(pos==i)) continue;end
%         v = res' * matrix(:, i);
%         if (i == 1)
%             mv = v;
%             p = i;
%         else
%             if (v > mv)
%                 mv = v;
%                 p = i;
%             end
%         end        
%     end
%     pos(c) = p;
%     res = mod(res + matrix(:, p), 2);
% end
% r(pos) = 1;