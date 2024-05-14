function unew = Lobatto_IIIC_4_Solver_B(u,A,dx,tau)

I = eye(size(A));
u(1)=0;u(end)=0;
u(2)=0;u(end-1)=0;
% u(3)=0;u(end-2)=0;
RHS = @(u) - A * u + computeF(u, dx);
RHS_prime = @(u) - A + computeF_prime(u,dx);

sqrt5 = sqrt(5);
a11 = 1/12;a12=-sqrt5/12;a13=sqrt5/12;a14=-1/12;
a21 = 1/12;a22=1/4;a23=(10-7*sqrt5)/60;a24=sqrt5/60;
a31 = 1/12;a32=(10+7*sqrt5)/60;a33=1/4;a34=-sqrt5/60;
a41 = 1/12;a42=5/12;a43=5/12;a44=1/12;

b1=1/12;b2=5/12;b3=5/12;b4=1/12;


        ques = @(uu) [ ...
            u + tau*a11*RHS(uu(1:end/4)) + tau*a12*RHS(uu((end/4+1):(2*end/4))) + tau*a13*RHS(uu((2*end/4+1):(3*end/4))) + tau*a14*RHS(uu((3*end/4+1):end)) - uu(1:end/4); ...
            u + tau*a21*RHS(uu(1:end/4)) + tau*a22*RHS(uu((end/4+1):(2*end/4))) + tau*a23*RHS(uu((2*end/4+1):(3*end/4))) + tau*a24*RHS(uu((3*end/4+1):end)) - uu((end/4+1):(2*end/4)); ...
            u + tau*a31*RHS(uu(1:end/4)) + tau*a32*RHS(uu((end/4+1):(2*end/4))) + tau*a33*RHS(uu((2*end/4+1):(3*end/4))) + tau*a34*RHS(uu((3*end/4+1):end)) - uu((2*end/4+1):(3*end/4)); ...
            u + tau*a41*RHS(uu(1:end/4)) + tau*a42*RHS(uu((end/4+1):(2*end/4))) + tau*a43*RHS(uu((2*end/4+1):(3*end/4))) + tau*a44*RHS(uu((3*end/4+1):end)) - uu((3*end/4+1):end)] ;
        ques_prime = @(uu) [ ...
            tau*a11*RHS_prime(uu(1:end/4)) - I, tau*a12*RHS_prime(uu((end/4+1):(2*end/4))), tau*a13*RHS_prime(uu((2*end/4+1):(3*end/4))), tau*a14*RHS_prime(uu((3*end/4+1):end));...
            tau*a21*RHS_prime(uu(1:end/4)), tau*a22*RHS_prime(uu((end/4+1):(2*end/4))) - I, tau*a23*RHS_prime(uu((2*end/4+1):(3*end/4))), tau*a24*RHS_prime(uu((3*end/4+1):end));...
            tau*a31*RHS_prime(uu(1:end/4)), tau*a32*RHS_prime(uu((end/4+1):(2*end/4))), tau*a33*RHS_prime(uu((2*end/4+1):(3*end/4))) - I, tau*a34*RHS_prime(uu((3*end/4+1):end)); ...
            tau*a41*RHS_prime(uu(1:end/4)), tau*a42*RHS_prime(uu((end/4+1):(2*end/4))), tau*a43*RHS_prime(uu((2*end/4+1):(3*end/4))), tau*a44*RHS_prime(uu((3*end/4+1):end)) - I];
        x0 = [u;u;u;u];
        answer = Newton_iter(ques, ques_prime, x0, 1e-12);

        u = u + tau*b1*RHS(answer(1:end/4)) + tau*b2*RHS(answer((end/4+1):(2*end/4))) + tau*b3*RHS(answer((2*end/4+1):(3*end/4))) + tau*b4*RHS(answer((3*end/4+1):end));

unew = u;
end
