% WCS: Select stage
% OUTPUT:
%          samplex --- 1D indices along hilbert curve
%          samplep --- 2D coordinates

function [samplex samplep] = wcs_select(cost, m, me, extra)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FOR DEBUG
VD = 0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flag2 = 0;
si = size(cost);
if (si(1) > 1 && si(2) > 1)
    flag2 = 1;
end
n = si(1) * si(2);

samplex = zeros(m, 1);
samplep = zeros(m, 2);

if (flag2 == 1) % 2D case
    costv = hilbertcurve(cost);
else 
    costv = cost;
end

if (strcmp(me, 'GREEDY'))
    [~, in] = sort(costv, 'ascend');
    samplex = in(1 : m);
    if (flag2 == 1)
        for i = 1 : m
            [samplep(i, 1) samplep(i, 2)] = invhcindex(samplex(i), si(1));
        end
    end
    
elseif (strcmp(me, 'C-A'))
    posb = extra * 1 ./ cost / sum(sum(1 ./ cost)) + (1 - extra) * ones(si(1), si(2)) / si(1) / si(2);
    
    sc = 0;
    i = 1;
    j = 1;
    bkposb = posb;
    while (sc < m)
        if ((bkposb(i, j) > 0) && (bkposb(i, j) >= 1 || binornd(1, bkposb(i, j)) == 1)) 
            sc = sc+1;
            samplep(sc, 1) = i;
            samplep(sc, 2) = j;
            samplex(sc) = hcindex(si(1), i, j);
            posb = posb / (1 - posb(i, j));
            posb(i, j) = 0;
            bkposb = posb;
        elseif (bkposb(i, j) > 0)
            bkposb = bkposb / (1 - bkposb(i, j));
            bkposb(i, j) = 0;
%             bkposb
%             pause;
            if (any(any(bkposb > 1)) == 1 || any(any(bkposb < 0)) == 1 && VD)
                bkposb
            end
        end
%         [xmax xin] = max(bkposb);
%         [ymax yin] = max(xmax);
%         j = yin;
%         i = xin(yin);
        j = j + 1;
        if (j > si(2))
            j = 1;
            i = i + 1;
        end
        if (i > si(1))
            i = 1;
        end
    end
    if (VD)
        figure
        plot(samplep(:,1), samplep(:, 2), 'o');
    end
    
elseif (strcmp(me, 'LLC'))
    xc = sqrt(m * si(1) / si(2));
    yc = sqrt(m * si(2) / si(1));
    xstep = si(1) / xc;
    ystep = si(2) / yc;
    sc = 0;
    
    if (floor(xc) < 1)
        xc = 1;
    end
    if (floor(yc) < 1)
        yc = 1;
    end

    for i = 1 : floor(xc)
        for j = 1 : floor(yc)
            sc = sc + 1;
            border = [floor((i - 1) * xstep + 1), floor(i * xstep), floor((j - 1) * ystep + 1), floor(j * ystep)];
            v = cost(border(1), border(3));
            samplep(sc, 1) = border(1);
            samplep(sc, 2) = border(3);
            samplex(sc) = hcindex(si(1), border(1), border(3));
            for xoff = border(1) : border(2)
                for yoff = border(3) : border(4)
                    if (cost(xoff, yoff) < v)
                        v = cost(xoff, yoff);
                        samplep(sc, 1) = xoff;
                        samplep(sc, 2) = yoff;
                        samplex(sc) = hcindex(si(1), xoff, yoff);
%                         fprintf('samplex(%i): %i\n', sc, samplex(sc));
                    end
                end
            end
        end
    end

    if (sc < m)
        randp = randperm(n);
        for i = 1 : n
            if (any(samplex == randp(i)) == 0)
                sc = sc + 1;
                samplex(sc) = randp(i);
                [samplep(sc, 1) samplep(sc, 2)] = invhcindex(randp(i), si(1));
                if (sc == m)
                    break;
                end
            end
        end
    end
    
elseif (strcmp(me, 'EVO'))
    info = ones(si(1), si(2));
    ipc = zeros (si(1), si(2));
    for i = 1 : si(1)
        for j = 1 : si(2)
            if (cost(i, j) ~= 0)
                ipc(i, j) = info(i, j) / cost(i, j);
            else
                ipc(i, j) = -1;
            end
        end
    end

    if (VD)
        ipcmonitor = figure;
        samplemonitor = figure;
        fademonitor = figure;
    end

    for sc = 1 : m
        if (any(ipc == -1) > 0)
            ma = -1;
        else
            ma = max(max(ipc));
        end

        candi = [];
        cn = 0;
        for i = 1 : si(1)
            for j = 1 : si(2)
                if (ipc(i, j) == ma)
                    cn = cn + 1;
                    candi(cn, 1) = i;
                    candi(cn, 2) = j;
                end
            end
        end

        elect = ceil(rand() * cn);
        samplep(sc, 1) = candi(elect, 1);
        samplep(sc, 2) = candi(elect, 2);
        samplex(sc) = hcindex(si(1), samplep(sc, 1), samplep(sc, 2));
    %     fprintf('Sample %d: (%d, %d) IPC=%.6f out of %d\n', sc, samplep(sc, 1), samplep(sc, 2), ipc(samplep(sc, 1), samplep(sc, 2)), length(candi));

        % multinomial
        sgm = 10;
        ampli = 5000;
        [X,Y] = meshgrid(1:si(1), 1:si(2));
        fade = ampli * mvnpdf([X(:), Y(:)], [samplep(sc, 1), samplep(sc, 2)], [sgm 0; 0 sgm]) + ones(n, 1);
        fade = reshape(fade, [si(1) si(2)]);

        ipc(samplep(sc, 1), samplep(sc, 2)) = 0;
        ipc = ipc ./ fade;

        if (VD)
            figure(samplemonitor);
            hold off
            contour(ipc);
            hold on
            plot(samplep(1:sc, 1), samplep(1:sc, 2), 'r^', 'MarkerFaceColor', 'r');
            
            figure(ipcmonitor);
            surf(ipc);
            
            figure(fademonitor);
            surf(fade);
            
%             pause(0.2);
        end
    end
end