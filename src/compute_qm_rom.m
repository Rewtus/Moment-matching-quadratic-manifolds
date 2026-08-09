function [V1, V2, W] = compute_qm_rom(A, B, S, L_1, L_2,r)
    % Compute V1 and V2 by solving Sylvester equations (QM Method)
    V1 = sylvester(-A, S, B*L_1);
    S_2 = (kron(S,eye(r))+kron(eye(r),S));
    V2 = sylvester(-A, S_2, B*L_2);
    W = V1; % Galerkin projection
    %V1,V2,W can be formed explicitly when S is diag
end