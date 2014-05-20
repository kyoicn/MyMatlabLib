%% simple interface to call KSVD

function [d, o] = simple_ksvd(data, N, sparsity, display)
if nargin < 4
    display = 0;
end
param.K = N;
param.L = sparsity;
param.numIteration = 50;
param.errorFlag = 0; % sparsity constraint
param.preserveDCAtom = 0;
% param.InitializationMethod = 'GivenMatrix';
% param.initialDictionary = rand(size(data, 1), N);
param.InitializationMethod = 'DataElements';
param.displayProgress = display;

[d, o] = KSVD(data, param);
end