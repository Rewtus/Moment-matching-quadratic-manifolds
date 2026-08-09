function [A, B, C] = wave_eqn(n, c, L, gamma)
    % Generates state-space matrices for the 1D Wave Equation
    %
    % Inputs:
    %   n - Number of spatial discretization points
    %   c - Wave speed 
    %   L - Length of the domain 
    %
    % Outputs:
    %   A, B, C - State-space matrices for dx/dt = Ax + Bu, y = Cx 

    % 1. Spatial Discretization
    dx = L / (n + 1);
    
    % 2. Construct the Discrete Laplacian (Finite Difference)
    main_diag = -2 * ones(n, 1);
    off_diag  = ones(n-1, 1);
    Laplacian = (1/dx^2) * (diag(main_diag) + diag(off_diag, 1) + diag(off_diag, -1));

    % 3. Assemble State-Space Matrices
    % State vector x = [displacement; velocity] -> size (2n x 1)
    
    % A = [ 0   I ]
    %     [ c^2*L  0 ]

    I = eye(n);
    Z = zeros(n);
    A = [zeros(n), eye(n); 
            (c^2 * Laplacian), -gamma * eye(n)];

    %(Assume force input at the center point)
    B = zeros(2*n, 1);
    B(round(n/2) + n) = 1; % Applying force to the 'velocity' state of middle node

    % Measure output as displacement at center point
    C = zeros(1,2*n);
    C(n/2)=1;

end