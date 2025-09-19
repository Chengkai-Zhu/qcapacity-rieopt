% This script computes the continuity upper bound in Theorem 1.
%
% Input: 
% rhoAB: The target state rhoAB.
% dA: Dimension of the system A.
% dB: Dimension of the system B.
% 
% Output:
% q_bound: The upper bound on D_->.
%
% This code is based on Theorem 1 in the paper.
% 
% (c) 2025, Chengkai Zhu.

function q_bound = DE_bound_isometry(rhoAB, dA, dB)

rhoAE = ComplementaryMap(rhoAB,[dA dB]); % the complementary state
s = size(rhoAE);
dE = s(1)/dA; % dimension of environment

cvx_begin sdp quiet
    variable J(dB*dE,dB*dE) hermitian;
    sigmaAE = PartialTrace(kron(eye(dA),J)*kron(PartialTranspose(rhoAB,2,[dA,dB]),eye(dE)),2,[dA,dB,dE]);
    eps = TraceNorm(sigmaAE - rhoAE)/2;
    minimize eps
    subject to
            J >= 0;
            PartialTrace(J,2,[dB,dE]) == eye(dB);
cvx_end

UBEG = chanconv(J, 'choi', 'isom', [dB dE]);
rhoAEG = kron(eye(dA), UBEG) * rhoAB * kron(eye(dA), UBEG');
sU = size(UBEG);
dG = sU(1)/dE;
rhoEG = PartialTrace(rhoAEG, 1, [dA sU(1)]);

HEG = Entropy(rhoEG) - Entropy(PartialTrace(rhoEG, 2, [dE dG]));

bi_entropy = - eps*log2(eps) - (1-eps)*log2((1-eps));

q_bound = HEG + eps*log2(dE^2-1)+ bi_entropy;  % bound via the epsilon degradable state
end

