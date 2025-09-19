% This function compute the cost function for the upper bound on D_->.
% 
% Input:
% X: The isometry that parameterizes the state extension.
% phiABE: The canonical purification of the target state.
% dA: Dimension of the system A.
% dB: Dimension of the system B.
% 
% Output:
% 
function f = bound_cost(X,phiABE,dA,dB)
rhoABFR = kron(eye(dA*dB), X) * phiABE * kron(eye(dA*dB), X'); % E to FR
sABFR = size(rhoABFR, 1);
dF = 2;
dR = sABFR/dA/dB/dF;
rhoABF = PartialTrace(rhoABFR, 2, [dA*dB*dF, dR]);
f = DE_bound_isometry(rhoABF);
end