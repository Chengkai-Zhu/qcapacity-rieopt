# qcapacity-rieopt

MATLAB code for the paper "+++"

This repository provides MATLAB code used to compute the upper and lower bounds on one-way distillable entanglement and quantum capacity based on the Riemannian optimization method developed in the paper ..


## Requirements

- QETLAB
- CVX
- Manopt
- cvxquad
- quantinf


## Files

### lowerbound

This folder includes the code for computing a lower bound on the quantum capacity via Riemannian gradient descent algorithm on a product of unitaris.

- `opt_localU.m`: The main file to do Riemannian gradient descent.

- `cohinfo_cost_local.m`: compute the cost function in terms of local unitaries.

- `cohinfo_grad_local.m`: compute the Euclidean gradient in terms of local unitaries.

- `gradient_wrapper.m`: help function to make the gradient struct.


### isotropic extension

This folder includes all extensions of qubit isotropic states we found that yield the best-known upper bound on the one-way distillable entanglement.
