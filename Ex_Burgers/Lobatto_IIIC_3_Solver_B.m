function unew = Lobatto_IIIC_3_Solver_B(u,A,dx,tau)

I = eye(size(A));
u(1)=0;u(end)=0;
u(2)=0;u(end-1)=0;
% u(3)=0;u(end-2)=0;
RHS = @(u) - A * u + computeF(u, dx);
RHS_prime = @(u) - A + computeF_prime(u,dx);


a11 = 1/6;a12=-1/3;a13=1/6;a21=1/6;a22=5/12;a23=-1/12;a31=1/6;a32=2/3;a33=1/6;
b1=1/6;b2=2/3;b3=1/6;

    
    ques = @(uu) [ ...
        u + tau*a11*RHS(uu(1:end/3)) + tau*a12*RHS(uu((end/3+1):(2*end/3))) + tau*a13*RHS(uu((2*end/3+1):end)) - uu(1:end/3); ...
        u + tau*a21*RHS(uu(1:end/3)) + tau*a22*RHS(uu((end/3+1):(2*end/3))) + tau*a23*RHS(uu((2*end/3+1):end)) - uu((end/3+1):(2*end/3)); ...
        u + tau*a31*RHS(uu(1:end/3)) + tau*a32*RHS(uu((end/3+1):(2*end/3))) + tau*a33*RHS(uu((2*end/3+1):end)) - uu((2*end/3+1):end)] ;
    ques_prime = @(uu) [ ...
        tau*a11*RHS_prime(uu(1:end/3)) - I, tau*a12*RHS_prime(uu((end/3+1):(2*end/3))), tau*a13*RHS_prime(uu((2*end/3+1):end));...
        tau*a21*RHS_prime(uu(1:end/3)), tau*a22*RHS_prime(uu((end/3+1):(2*end/3))) - I, tau*a23*RHS_prime(uu((2*end/3+1):end));...
        tau*a31*RHS_prime(uu(1:end/3)), tau*a32*RHS_prime(uu((end/3+1):(2*end/3))), tau*a33*RHS_prime(uu((2*end/3+1):end)) - I];
    x0 = [u;u;u];
    answer = Newton_iter(ques, ques_prime, x0, 1e-12);
    
    u = u + tau*b1*RHS(answer(1:end/3)) + tau*b2*RHS(answer((end/3+1):(2*end/3))) + tau*b3*RHS(answer((2*end/3+1):end));

unew = u;
end
