function rho_grouped = ncopy_interleave(rhoAB, dA, dB, n)

    rho_total = 1;
    for i = 1:n
        rho_total = kron(rho_total, rhoAB);
    end

    if n == 1
        rho_grouped = rho_total;
        return;
    end

    num_systems = 2 * n;
    p_odd = 1:2:num_systems;  
    p_even = 2:2:num_systems; 
    perm = [p_odd, p_even]; 

    dims = repmat([dA, dB], 1, n);

    rho_grouped = PermuteSystems(rho_total, perm, dims);
end


function CI = CohInfo(state,da,db)
rho_BF = PartialTrace(state, 1, [da,db]);
CI = Entropy(rho_BF) - Entropy(state); % coherant information
end