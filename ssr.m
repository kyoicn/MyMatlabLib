%% Singal to signal ratio
% s1: 
% s2: 

function r=ssr(s1,s2)
si1 = size(s1);
si2 = size(s2);

if (si1(1) == si2(2) && si1(2) == si2(1)) %transpose
    s2 = s2';
elseif (si1(1) ~= si2(1) || si1(2) ~= si2(2))
    fprintf('Signal dimensions must agree.\n');
    r = 0;
    return;
end
r=10*log10(vnorm(s1,2)^2/vnorm(s2,2)^2);