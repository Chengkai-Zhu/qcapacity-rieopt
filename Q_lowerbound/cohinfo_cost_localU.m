% The function for calculating the cost function of channel coherent information
%
% Input: 
% X: A product of unitaries.
% n: Number of copies of the channel.
% d: Local dimension of a bipartite state (same dim).
% dR: Dimension of the auxiliary system.
% K: A cell array containing the Kraus operators {K_i} for the n-copy channel.
%
% Output: The cost function.
% 
% This code is based on Algorithm 3 in the paper.
% 
% (c) 2025, Chengkai Zhu.

function f = cohinfo_cost_localU(X, n, d, dR, K)
% calculate the coherent information
num_U = 2*n - 1;
dims = [dR, repmat(d, 1, n)]; % Vector of dimensions for each system
total_dim = prod(dims);

% Initial state |0>_{RA^n}
p0 = zeros(total_dim, 1);
p0(1) = 1;
phi = kron(X.R1, eye(total_dim/dR/d)) * p0;

for k = 2:num_U-1
    U_k = X.Main(:,:,k-1);
    % Determine which two systems U_k acts on
    if k <= n
        sys_idx = [k, k+1];
    else
        j = 2*n - k;
        sys_idx = [j, j+1];
    end
    % Apply the action of V_k = U_k \otimes I on phi{k}
    phi = apply_local_unitary(phi, U_k, sys_idx, dims);
end
phi = kron(X.R2, eye(total_dim/dR/d)) * phi;

rho = 0;
for i=1:length(K)
    IK = kron(eye(dR), K{i});
    psi = IK*phi;
    rho = rho + psi*psi';
end

s = size(K{1});
dout = s(1);
rho_B = PartialTrace(rho, 1, [dR,dout]);
cohinfo = Entropy(rho_B) - Entropy(rho);
f = -cohinfo;
end