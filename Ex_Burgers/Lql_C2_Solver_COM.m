function unew = Lql_C2_Solver_COM(u, A, dx,dt)
    u(1)=0;u(end)=0;
    u(2)=0;u(end-1)=0;
    A=dt*A;
    I = eye(size(A));
    f1 = computeF(u, dx);
       
    I_Euler_op = (3.149*A^2 + 0.775*A + I);      
    [L_IE,U_IE] = lu(I_Euler_op);
    
     r = (3.125*A + I) * f1;
     b = (I-0.225*A + 0.024*A^2)* u + dt * r;
     unew = (U_IE\(L_IE\(b)));   
end
