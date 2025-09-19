% The function for canonical purification of a state
%
% Input: Target state rhoAB.
% Output: 
% phiABE: The purification phi_ABE.
% r: Dimension of the auxiliary system.

function [phi_ABE, r] = canoPurif(rho_AB)
dAB = size(rho_AB, 1);
rho_AB = (rho_AB + rho_AB') / 2;

% Diagonalize rho_AB to get eigenvalues and eigenvectors
[V, D] = eig(rho_AB, 'vector'); % 'vector' returns eigenvalues as a vector
[p, idx] = sort(D, 'descend');   % Sort eigenvalues in descending order
V = V(:, idx);                   % Sort eigenvectors accordingly

% Identify non-zero eigenvalues above tolerance
tol = 1e-12;
non_zero = p > tol;
r = sum(non_zero);               % Rank of rho_AB (dimension of E)
p_nonzero = p(non_zero);         % Non-zero eigenvalues
V_nonzero = V(:, non_zero);      % Corresponding eigenvectors

% Initialize the purified state vector
phi_ABE = zeros(dAB * r, 1);

for i = 1:r
    psi_AB = V_nonzero(:, i);
    % Basis vector for environment E (standard basis in r dimensions)
    e_i = zeros(r, 1);
    e_i(i) = 1;
    % Tensor product: psi_AB (AB) x e_i (E)
    phi_ABE = phi_ABE + sqrt(p_nonzero(i)) * kron(psi_AB, e_i);
end
end