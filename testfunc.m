%OMP Recovery
%
% stop: stop condition

function re = testfunc(a)
re = a + 1;
a = 0;
re = re + a;