function [hsv] = compute_hsv(A, B, C)
   
    % Compute Cholesky factors of Gramians
    U = lyapchol(A, B);   
    L = lyapchol(A', C'); 
    
    % SVD of the product of Cholesky factors
    [~, Sigma, ~] = svd(L * U');
    
    % Output the Hankel Singular Values as a vector
    hsv = diag(Sigma);
end