function unew = Lobatto_IIIC_2_Solver_B(u,A,dx,tau)

I = eye(size(A));
u(1)=0;u(end)=0;
u(2)=0;u(end-1)=0;
% u(3)=0;u(end-2)=0;
RHS = @(u) - A * u + computeF(u, dx);
RHS_prime = @(u) - A + computeF_prime(u,dx);


a11 = 1/2;a12 = -1/2;a21 = 1/2;a22 = 1/2;
b1 = 1/2; b2 = 1/2;

    
    ques = @(uu) [u + tau*a11*RHS(uu(1:end/2)) + tau*a12*RHS(uu((end/2+1):end)) - uu(1:end/2); ...
        u + tau*a21*RHS(uu(1:end/2)) + tau*a22*RHS(uu((end/2+1):end)) - uu((end/2+1):end)] ;
    ques_prime = @(uu) [tau*a11*RHS_prime(uu(1:end/2)) - I, tau*a12*RHS_prime(uu((end/2+1):end));
        tau*a21*RHS_prime(uu(1:end/2)), tau*a22*RHS_prime(uu((end/2+1):end)) - I];
    x0 = [u;u];
    answer = Newton_iter(ques, ques_prime, x0, 1e-12);
    
    u = u + tau*b1*RHS(answer(1:end/2)) + tau*b2*RHS(answer((end/2+1):end));

   
unew = u;

end
