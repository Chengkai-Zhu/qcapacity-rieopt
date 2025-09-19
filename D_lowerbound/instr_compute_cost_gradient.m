% Computes the cost function in optimizing D^{(1)} over unitary manifold
% that parameterizes Kraus operators.
% System ordering convention: The full Hilbert space is ordered A, then B, then M.
% Input:
% X: The unitary for instrument.
% rho_ab: The target bipartite state.
% dimA: Dimension of the system A.
% dimB: Dimension of the system B.
% dimM: Dimension of the classical system M.
% 
% Output:
% cost: the coherent information of the output state.
% grad_U: the Euclidean gradient w.r.t. the unitary.
%
% This code is based on Algorithm 2 in the paper.
% 
% (c) 2025, Chengkai Zhu.

function [cost, grad_U] = instr_compute_cost_gradient(X, rho_ab, dimA, dimB, dimM)
    
    dimAB = dimA * dimB;
    dimABM = dimA * dimB * dimM;

    idA = eye(dimA);
    idB = eye(dimB);
    idAB = eye(dimAB);
    
    
    % Initial state on M is |0><0|
    ket0_M = [1; zeros(dimM-1, 1)];
    proj0_M = ket0_M * ket0_M';
    
    % Full initial state on ABM: rho_ab kron |0><0|_M
    rho0 = kron(rho_ab, proj0_M);
    % The unitary U_am acts on systems A and M. In our ABM ordering, these
    U_am_full = PermuteSystems(kron(X, idB), [1 3 2], [dimA, dimM, dimB]);

    % the state after the channel, Psi
    Psi = U_am_full * rho0 * U_am_full';
        
    sigma_ab_blocks = cell(1, dimM);
    sigma_b_blocks = cell(1, dimM);
    log_sigma_ab_blocks = cell(1, dimM);
    log_sigma_b_blocks = cell(1, dimM);
    
    cost_S_ab = 0;
    cost_S_b = 0;
    
    for j = 1:dimM
        % Projector for the j-th outcome on system M
        ket_j_M = zeros(dimM, 1);
        ket_j_M(j) = 1;
        
        % Operator to project M onto <j| and get a state on AB
        % This is equivalent to Tr_M(Psi * (Id_AB kron |j><j|_M))
        proj_op = kron(idAB, ket_j_M');
        sigma_ab_j = proj_op * Psi * proj_op';
        sigma_ab_blocks{j} = sigma_ab_j;
        
        % Compute sigma_b by tracing out system A
        sigma_b_j = PartialTrace(sigma_ab_j, 1, [dimA, dimB]);
        sigma_b_blocks{j} = sigma_b_j;
        
        % Add a small identity for numerical stability before taking log
        % stable_sigma_ab_j = sigma_ab_j + epsilon * eye(size(sigma_ab_j));
        % stable_sigma_b_j = sigma_b_j + epsilon * eye(size(sigma_b_j));
        
        log_sigma_ab_blocks{j} = logm(sigma_ab_j)/log(2);
        log_sigma_b_blocks{j} = logm(sigma_b_j)/log(2);
        
        % S(X) = Tr(X log X)
        cost_S_ab = cost_S_ab + trace(sigma_ab_j * log_sigma_ab_blocks{j});
        cost_S_b = cost_S_b + trace(sigma_b_j * log_sigma_b_blocks{j});
    end
    
    % cost function
    cost = real(cost_S_b - cost_S_ab);
    
    % construct the Y operators for the gradient calculation
    Y_2 = zeros(dimABM);
    Y_1 = zeros(dimABM);
    
    for j = 1:dimM
        ket_j_M = zeros(dimM, 1);
        ket_j_M(j) = 1;
        proj_j_M = ket_j_M * ket_j_M';
        
        % Y_2 is a block-diagonal operator with log(sigma_ab_j) on the j-th block
        Y_2 = Y_2 + kron(log_sigma_ab_blocks{j}, proj_j_M);
        
        % Y_1 is block-diagonal with (Id_A kron log(sigma_b_j)) on the j-th block
        Y_1 = Y_1 + kron(kron(idA, log_sigma_b_blocks{j}), proj_j_M);
    end
    
    Y = Y_1 - Y_2;
        
    % Core of the gradient formula: G = U_full' * Y * U_full * rho0
    G = Y * U_am_full * rho0;

    % Trace out system B (the 2nd system in our ABM order) to get a matrix on AM
    G_am = PartialTrace(G, 2, [dimA, dimB, dimM]);
    
    % The Euclidean gradient is the transpose of this result
    grad_U = 2*G_am;
end