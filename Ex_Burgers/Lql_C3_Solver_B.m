function unew = Lql_C3_Solver_B(u, A, dx,dt)
    u(1)=0;u(end)=0;
    u(2)=0;u(end-1)=0;
    A=dt*A;
    I = eye(size(A));
    f1 = computeF(u, dx);
    
    I_Euler_op = I+0.78986*A+0.38283*A*A;   
    [L_IE,U_IE] = lu(I_Euler_op);

    
    r = (0.37797*A+I) * f1;
    b = (I-0.21014*A + 0.00486*A*A)*u + dt * r;
    unew = (U_IE\(L_IE\b));    
    
  
    

end
