function uNew = Lql_C2I_Solver_B(u, A, dx, dt)
    maxIter = 40;
    tolerance = 1e-12;
    Nx = length(u) - 1;
    I = eye(Nx+1);
    A = dt * A;
    u(1)=0;u(end)=0;
    u(2)=0;u(end-1)=0;

%     I_Euler_op = I+A;
%     [L_IE, U_IE] = lu(I_Euler_op);
%     
%     % Calculate R and P functions
%     R = (U_IE \ (L_IE \ (I )));
%     P = (U_IE \ (L_IE \ (I )));
    
    I_Euler_op = I + 0.79033 * A + 0.37931 * A^2;
    [L_IE, U_IE] = lu(I_Euler_op);
    
    % Calculate R and P functions
    R = (U_IE \ (L_IE \ (I - 0.20967 * A + 0.00484 * A^2)));
    P = (U_IE \ (L_IE \ (I + 0.37447 * A)));
    
    % Initialization
    uNew = u;
    F = zeros(Nx+1,1);
    J = zeros(Nx+1,Nx+1);
    
     for iter = 1:maxIter
        unew(1)=0;unew(end)=0;
        % Evaluate F using loops
        for i = 2:Nx
            sum_Ru = 0;
            sum_Pu_diff = 0;
            for j = 2:Nx
                sum_Ru = sum_Ru + R(i,j) * u(j);
                sum_Pu_diff = sum_Pu_diff + P(i,j) * (uNew(j-1)^2 - uNew(j+1)^2);
            end
            F(i) = uNew(i) - sum_Ru - (dt / (4 * dx)) * sum_Pu_diff;
        end
        
        % Boundary conditions
        F(1) = 0;
        F(end) = 0;
        
        % Compute Jacobian using loops
        for i = 2:Nx
            for k = 2:Nx
                if k == i-1
                    J(i,k) =  - (dt / (2 * dx)) * P(i,k) * uNew(i-1);
                elseif k == i
                    J(i,k) = 1 ;
                elseif k == i+1
                    J(i,k) =  (dt / (2 * dx)) * P(i,k) * uNew(i+1);
                else
                    J(i,k) = 0;
                end
            end
        end
        
        % Boundary conditions for Jacobian
       J(1, :) = 0; J(end, :) = 0;
        J(1, 1) = 1; J(end, end) = 1;
        
        % Newton-Raphson update
        deltaU = J \ (-F);
        uNew = uNew + deltaU;
        
        % Check convergence
        if max(abs(deltaU)) < tolerance
            break;
        end
    end
end
