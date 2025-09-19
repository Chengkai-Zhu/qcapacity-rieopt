% The main demo script for optimizing the upper bound on D_-> of the isotropic states.

clear;
d=2;
p=0.1425/(3/4); % depolarizing parameter
bell_state = MaxEntangled(d)*MaxEntangled(d)';
id = eye(d^2)/d^2;
iso_state = (1-p)*bell_state + p*id;

[phiABE, dE]= canoPurif(iso_state);
phiABE = phiABE*phiABE';

dF = 2;

%% define manifold
problem.M = stiefelcomplexfactory(dE*dF,dE);
problem.cost = @(X) bound_cost(X,phiABE,d);

% check whether the gradient is correct
% checkgradient(problem)

options.maxiter = 1000;
options.tolgradnorm = 1e-10;

% choose optimization method
[Xopt, f, info] = steepestdescent(problem, Xopt, options);


