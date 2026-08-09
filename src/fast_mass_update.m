function M_update = fast_mass_update(E2r, x, r)
    % Evaluates E2r * (kron(eye(r), x) + kron(x, eye(r))) in strict O(r^3) time
    M_update = zeros(r, r);
    for i = 1:r
        % Extract the i-th block of size r x r from E2r
        E_i = E2r(:, (i-1)*r + 1 : i*r); 
        
        % 1. Action of (x \otimes I): scales and sums the blocks
        M_update = M_update + x(i) * E_i; 
        
        % 2. Action of (I \otimes x): forms the columns of the output
        M_update(:, i) = M_update(:, i) + E_i * x; 
    end
end