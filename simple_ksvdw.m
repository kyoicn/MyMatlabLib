%% simple interface to call KSVDW

function [d, o] = simple_ksvdw(data, N, sparsity, weights, display)
if nargin < 5
    display = 0;
elseif nargin < 4
    weights = ones(size(data,1), 1);
end
param.K = N;
param.L = sparsity;
param.numIteration = 100;
param.errorFlag = 0; % sparsity constraint
param.preserveDCAtom = 0;
param.InitializationMethod = 'DataElements';
param.displayProgress = display;

[d, o] = KSVDW(data, param, weights);
end