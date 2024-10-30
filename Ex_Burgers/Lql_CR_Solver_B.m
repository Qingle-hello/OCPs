function unew = Lql_CR_Solver_B(u, A, dx,dt)

    u(1)=0;u(end)=0;
    u(2)=0;u(end-1)=0;
    % u(3)=0;u(end-2)=0;

    A=dt*A;
    I = eye(size(A));
    f1 = computeF(u, dx);
    
    I_Euler_op = (0.37922*A^2 + 0.78875*A + I);      

    B = inv(I_Euler_op);
    
    r = (0.37443*A + I) * f1;
    b = (0.98737*I - 0.22121*A)* u + dt * r;
    unew = B * b + 0.01263*u;

  
end
