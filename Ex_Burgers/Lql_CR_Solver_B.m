function unew = Lql_CR_Solver_B(u, A, dx,dt)

    u(1)=0;u(end)=0;
    u(2)=0;u(end-1)=0;
    % u(3)=0;u(end-2)=0;

    A=dt*A;
    I = eye(size(A));
    f1 = computeF(u, dx);
    
    I_Euler_op = (0.37922*A^2 + 0.78875*A + I);      
    [L_IE,U_IE] = lu(I_Euler_op);
    
    r = (0.37443*A + I) * f1;
    b = (I - 0.21125*A + 0.00479*A^2)* u + dt * r;
    unew = (U_IE\(L_IE\(b)));   

  
end
