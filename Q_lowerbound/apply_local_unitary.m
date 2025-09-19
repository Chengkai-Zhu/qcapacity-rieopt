% Applies a local two-system unitary U to a state vector phi_in.
%
% Input:
% phi_in: Input state vector.
% U: The (d*d)x(d*d) unitary matrix.
% target_sys: A 1x2 vector with the indices of the systems U acts on.
% dims: A vector of dimensions for all subsystems.
%
% Output: the code state.
%
% This code is based on Algorithm 3 in the paper.
% 
% (c) 2025, Chengkai Zhu.

function phi_out = apply_local_unitary(phi_in, U, target_sys, dims)

    num_systems = length(dims);
    total_dim = prod(dims);
    
    % Permute target systems to the front
    other_sys = setdiff(1:num_systems, target_sys);
    perm_fwd = [target_sys, other_sys];
    phi_perm = PermuteSystems(phi_in, perm_fwd, dims);
    
    % apply unitary
    d_target = prod(dims(target_sys));
    d_other = total_dim / d_target;
    
    phi_out = kron(U, eye(d_other)) * phi_perm;
    
    % permute systems to original order
    [~, perm_bwd] = sort(perm_fwd); % Inverse permutation
    
    phi_out = PermuteSystems(phi_out, perm_bwd, dims(perm_fwd));
end