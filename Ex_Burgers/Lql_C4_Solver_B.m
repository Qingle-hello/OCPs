function unew = Lql_C4_Solver_B(u, A, dx,dt)

    u(1)=0;u(end)=0;
    u(2)=0;u(end-1)=0;
    A=dt*A;
    I = eye(size(A));
    f1 = computeF(u, dx);
    
    I_Euler_op = (0.37898*A^2 + 0.78909*A + I);  
    B = inv(I_Euler_op);

    
    r = (0.37422*A + I) * f1;
    b = (0.98744*I - 0.22082*A)* u + dt * r;

    unew = B*b + 0.01256 * u;
    
end
