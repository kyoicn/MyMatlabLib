%Extended MCR

function re = exmcr(s, k, ex)
n=length(s);
s=s';
%with extended
for i = n + 1 : n + ex
    s = [s s(n)];
end
[t, b] = dct4(s);
[a, phi] = qr(b(1 : k, 1 : (n + ex)));
a = a';
sink_gets = phi * s';
sol = kramer(a(1 : k, 1 : k), sink_gets(1 : k));
for i = 1 : n + ex
    if i <= k
        exr(i) = sol(i);
    else
        exr(i) = 0;
    end
end
re = b * exr';
re=re(1:n);