# qcapacity-rieopt

MATLAB code for the paper "Geometric optimization for quantum communication"

This repository provides MATLAB code used to compute the upper and lower bounds on one-way distillable entanglement and quantum capacity based on the Riemannian optimization method developed in the paper ..


## Requirements

- QETLAB
- CVX
- Manopt
- cvxquad
- quantinf


## Files

### D lowerbound

This folder includes the code for computing a lower bound on the one-way distillable entanglement via Riemannian gradient descent algorithm on a unitary manifold.

- `opt_instr.m`: The main function to do Riemannian gradient descent.

- `instr_compute_cost_gradient.m`: compute the cost function and the gradient in terms of the local unitary that parameterize the instrument.

- `ncopy_interleave.m`: prepare the n-copy bipartite state with A^n|B^n partition.


### Q lowerbound

This folder includes the code for computing a lower bound on the quantum capacity via Riemannian gradient descent algorithm on a product of unitary manifolds.

- `opt_localU.m`: The main function to do Riemannian gradient descent.

- `cohinfo_cost_local.m`: compute the cost function in terms of local unitaries.

- `cohinfo_grad_local.m`: compute the Euclidean gradient in terms of local unitaries.

- `gradient_wrapper.m`: help function to make the gradient struct.

### DamDeph codestates

This folder includes all optimized code states for the damping-dephasing channel.

### GADC codestates

This folder includes all optimized code states for the generalized amplitude damping channel.

### isotropic extensions

This folder includes all extensions of qubit isotropic states we found that yield the best-known upper bound on the one-way distillable entanglement.
