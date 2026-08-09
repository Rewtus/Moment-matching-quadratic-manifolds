function [A, B, C] = cde(n, L, v)
% Generates a 1D LTI system for linear transport eqn
%
%   Inputs:
%       n - Number of states (interior + outlet nodes)
%       L - Length of the spatial domain
%       v - velocity chosen positive (right)
%
%   Outputs:
%       A - System matrix (n x n)
%       B - Input matrix (n x 1)
%       C - Output matrix (1 x n)

    % Total nodes including the right boundary element z=L
    N = n; 
    
    % Node 1 is z=1 (outlet), Node N+1 is z=-1 (inlet)
    x = cos(pi * (0:N)' / N);
    
    c = [2; ones(N-1, 1); 2];

    D = zeros(N+1, N+1);
    for r = 1:N+1
        for col = 1:N+1
            if r ~= col
                D(r,col) = (c(r)/c(col)) * ((-1)^((r-1)+(col-1))) / (x(r) - x(col));
            end
        end
    end

    for r = 1:N+1
        D(r,r) = -sum(D(r,:));
    end
    
    D_physical = (2 / L) * D;
    
    % Node N+1 is upstream inlet boundary condition: x(N+1) = u(t)

    A = -v * D_physical(1:N, 1:N);
    B = -v * D_physical(1:N, N+1);
    
    % Measure at the physical exit z = L 
    C = zeros(1, n);
    C(1) = 1;
end