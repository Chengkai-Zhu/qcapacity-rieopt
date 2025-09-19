% The main script for optimizing the unitaries.
%
% Input:
% n: Number of copies of the channel
% dR: dimension of the auxiliary system
% K: a cell array containing the Kraus operators {K_i} for the n-copy channel.
%
% Output:
% Xopt: optimized unitaries
% optbound: lower bound on Q^{(1)}(N^n)
%
% This code is based on Algorithm 3 in the paper.
% 
% (c) 2025, Chengkai Zhu.

function [Xopt, optbound] = opt_localU(n, dR, K)

d = 2;
% n = 6; % n copies of channels
% dR = 2;

s0 = [1 0;0 0];
s1 = [0 0;0 1];

K_ncopy = NKraus(K, n); % Kraus operators of the n-copy channel


%% define manifold
num_unitary = 2*n-1;

elements.R1 = unitaryfactory(dR*d, 1);

elements.Main = unitaryfactory(d^2, num_unitary-2);

elements.R2 = unitaryfactory(dR*d, 1);

problem.M = productmanifold(elements);

problem.cost = @(X) cohinfo_cost_localU(X, n, d, dR, K_ncopy);

problem.egrad = @(X) gradient_wrapper(X, n, d, dR, K_ncopy);

% check whether the gradient is correct
% checkgradient(problem)

options.maxiter = 1000;
options.tolgradnorm = 1e-7;
options.verbosity = 0;

% choose optimization method
[Xopt, f, ~] = conjugategradient(problem, [], options);


optbound = -f/n;
% fprintf('max coh info: %f\n', -f/n);
end









