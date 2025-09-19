% The main script for optimizing the instrument.
%
% Inputs:
% rhoAB: The target state.
% n: Number of copies of the state.
% da: Dimension of the system A.
% da: Dimension of the system B.
% dm: Dimension of the classical system M.
%
% Output:
% Xopt: optimimal unitary for the instrument.
% optbound: lower bound on D^{(1)}(rho^n).
%
% This code is based on Algorithm 2 in the paper.
% 
% (c) 2025, Chengkai Zhu.

function [Xopt, optbound] = opt_instr(rhoAB, da, db, dm, n)

% define Stiefel manifold
da = d;
db = d;
dm = 2;
n = 5;

rhoAB_ncopy = n_copy_interleave(rhoAB, da, db, n);

%% define manifold
problem.M = unitaryfactory(da^n*dm);
problem.cost = @(X) get_cost(X, rhoAB_ncopy, da^n, db^n, dm);
problem.egrad = @(X) get_grad(X, rhoAB_ncopy, da^n, db^n, dm);

% check whether the gradient is correct
% checkgradient(problem)xia

options.maxiter = 200;
options.tolgradnorm = 1e-7;

% choose optimization method
% [Xopt, f, info] = steepestdescent(problem, [], options);
[Xopt, f, info] = conjugategradient(problem, [], options);

optbound = -f/n;
fprintf('max coh info: %f\n', -f/n);

end


function cost = get_cost(X, rho_ab, dimA, dimB, dimM)
[cost, ~] = instr_compute_cost_gradient(X, rho_ab, dimA, dimB, dimM);
end

function grad_U = get_grad(X, rho_ab, dimA, dimB, dimM)
[~, grad_U] = instr_compute_cost_gradient(X, rho_ab, dimA, dimB, dimM);
end
