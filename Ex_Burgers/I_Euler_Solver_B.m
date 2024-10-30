function unew = I_Euler_Solver_B(u, A, dx,dt)

    u(1)=0;u(end)=0;
    u(2)=0;u(end-1)=0;
    A=dt*A;
    I = eye(size(A));
    f1 = computeF(u, dx);
    
    I_Euler_op = (I+A);    
    B = inv(I_Euler_op);


    unew = B*(u + dt * f1 );     

    

end
