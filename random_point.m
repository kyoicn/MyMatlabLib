%% pick up random points
% n the number of points
% d1 dimension 1
% d2 dimension 2

function re=random_point(n, d1, d2)
re = zeros(2, n);
c = 0;
for i=1:1:n
    flag = 0;
    while (flag == 0)
        x=unidrnd(d1);
        y=unidrnd(d2);
        flag = 1;
        for j = 1:1:c
            if (re(1,j) == x & re(2,j) == y)
                flag = 0;
                break;
            end
        end
        re(1,i) = x;
        re(2,i) = y;
    end
end