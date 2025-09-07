function [Xopt, optbound] = GAD_localU(n, dR)

d = 2;
% n = 6; % n copies of channels
% dR = 2;

s0 = [1 0;0 0];
s1 = [0 0;0 1];

% generalized amplitude damping channel
gamma = 0.44035;
N = 0.1;
KAD{1} = sqrt(1-N)*(s0 + sqrt(1-gamma)*s1);
KAD{2} = sqrt(gamma*(1-N))*[0 1;0 0];
KAD{3} = sqrt(N)*(sqrt(1-gamma)*s0 + s1);
KAD{4} = sqrt(gamma*N)*[0 0;1 0];
KADncopy = NKraus(KAD, n); % Kraus operators of the n-copy depolarizing channel


%% define manifold
num_unitary = 2*n-1;

elements.R1 = unitaryfactory(dR*d, 1);

elements.Main = unitaryfactory(d^2, num_unitary-2);

elements.R2 = unitaryfactory(dR*d, 1);

problem.M = productmanifold(elements);

problem.cost = @(X) cohinfo_cost_localU(X, n, d, dR, KADncopy);

problem.egrad = @(X) gradient_wrapper(X, n, d, dR, KADncopy);

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









