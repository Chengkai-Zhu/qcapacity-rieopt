% Computes the Euclidean gradient of the coherent information for a MERA-like circuit.
%
% Inputs:
% X: A (d^2) x (d^2) x (2n-1) array where X(:,:,k) is the k-th local unitary U_k.
% n: Number of systems in the A part of the chain.
% d: Dimension of each system in A (e.g., d=2 for qubits).
% dR: Dimension of the auxiliary system R.
% K: A cell array containing the Kraus operators {K_i} for the n-copy channel.
%
% Output:
% grad: A (d^2) x (d^2) x (2n-1) array of Euclidean gradients for each U_k.
%
% This code is based on Algorithm 3 in the paper.
% 
% (c) 2025, Chengkai Zhu.

function [gradR1, gradR2, grad] = cohinfo_grad_localU(X, n, d, dR, K)

num_U = 2*n - 1;
num_systems = n + 1; % Total systems: R, A_1, ..., A_n
dims = [dR, repmat(d, 1, n)]; % Vector of dimensions for each system
total_dim = prod(dims);

% Store intermediate states
phi = cell(1, num_U + 1);

% Initial state |0>_{RA^n}
p0 = zeros(total_dim, 1);
p0(1) = 1;
phi{1} = p0;
phi{2} = kron(X.R1, eye(total_dim/dR/d)) * phi{1};

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
    phi{k+1} = apply_local_unitary(phi{k}, U_k, sys_idx, dims);
end
phi{num_U+1} = kron(X.R2, eye(total_dim/dR/d)) * phi{num_U};

psi = phi{end};

% Compute the Gradient Signal 'G'
Npsi = 0;
for i=1:length(K)
    Npsi = Npsi + (kron(eye(dR), K{i}) * psi) * (psi' * kron(eye(dR), K{i}'));
end

s = size(K{1});
dout = s(1); % Output dimension of the channel's action space

% The remaining space has dimension dout = d^n.
NpsiB = logm(PartialTrace(Npsi, 1, [dR dout]))/log(2);

NNpsiB = 0;
for i=1:length(K)
    NNpsiB = NNpsiB + K{i}' * NpsiB * K{i};
end

logNpsi = logm(Npsi)/log(2);

NNpsi = 0;
for i=1:length(K)
    NNpsi = NNpsi + kron(eye(dR),K{i}') * logNpsi * kron(eye(dR), K{i});
end

G = 2 * (kron(eye(dR), NNpsiB) - NNpsi) * psi;

% Gradient Calculation and Propagation
grad = zeros(d^2, d^2, num_U);

gradR2 = PartialTrace(G*phi{num_U}', 3:num_systems, dims);

G_back = apply_local_unitary(G, X.R2', [1 2], dims);

for k = num_U-1:-1:2
    U_k = X.Main(:,:,k-1);
    right_vec = phi{k}; % This is |\phi_{k-1}>
    
    % Determine which two systems U_k acted on
    if k <= n
        sys_idx = [k, k+1];
    else
        j = 2*n - k;
        sys_idx = [j, j+1];
    end
    
    grad(:,:,k) = PartialTrace(G_back*right_vec', setdiff(1:num_systems, sys_idx), dims);
    
    % This computes G_back_new = V_k' * G_back_old
    G_back = apply_local_unitary(G_back, U_k', sys_idx, dims);
end
    gradR1 = PartialTrace(G_back*phi{1}', 3:num_systems, dims);

end