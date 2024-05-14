function unew = Radau_IIA_3_Solver_B(u,A,dx,tau)

I = eye(size(A));
u(1)=0;u(end)=0;
u(2)=0;u(end-1)=0;
u(3)=0;u(end-2)=0;

RHS = @(u) - A * u + computeF(u, dx);
RHS_prime = @(u) - A + computeF_prime(u,dx);

s=sqrt(6);

a11 = (88-7*s)/360;a12=(296-169*s)/1800;a13=(-2+3*s)/225;
a21=(296+169*s)/1800;a22=(88+7*s)/360;a23=(-2-3*s)/225;
a31=(16-s)/36;a32=(16+s)/36;a33=1/9;

b1=(16-s)/36;b2=(16+s)/36;b3=1/9;


    
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
