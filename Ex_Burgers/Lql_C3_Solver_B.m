function unew = Lql_C3_Solver_B(u, A, dx,dt)
    u(1)=0;u(end)=0;
    u(2)=0;u(end-1)=0;
    A=dt*A;
    I = eye(size(A));
    f1 = computeF(u, dx);
    
    I_Euler_op = I+0.78986*A+0.38283*A*A;   
    B = inv(I_Euler_op);
    
    r = (0.37797*A+I) * f1;
    b = (0.98730*I - 0.22017*A)*u + dt * r;
    unew = B*b + 0.0127 * u;
    
  
    

end
