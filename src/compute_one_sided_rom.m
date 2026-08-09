function [A_osr_r, B_osr_r, C_osr_r, V_orth] = compute_one_sided_rom(A, B, C, freqs, r)

    n=size(A,1);
    V_osr=zeros(n,r/2);

    %Explicitly form
    for i=1:r/2
        V_osr(:,i)=(freqs(i)*eye(n)-A)\B;
    end

    %Using conjugate pairs and enforcing realness
    V_real = [real(V_osr), imag(V_osr)];

    %Orthogonalize to prevent unstable A_r 
    %This results in smaller state space dimension 
    V_orth=orth(V_real); %Time taking step!

    A_osr_r = V_orth' * A * V_orth;
    B_osr_r = V_orth' * B;
    C_osr_r = C * V_orth;
end
