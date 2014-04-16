%%kramer solver

function r = kramer(a, b)
n = length(b);
size_a = size(a);
if (size_a(1) ~= n || size_a(2) ~= n)
    r = 0;
else
    r = [];
    basis = det(a);
    for i = 1 : n
        diff = [];
        for j = 1 : (i - 1)
            diff = [diff a(1 : n, j)];
        end
        diff = [diff b];
        for j = (i + 1) : n
            diff = [diff a(1 : n, j)];
        end
        r(i) = det(diff) / basis;
        diff = [];
    end
end