function [uend, time1, cut, energy] = Lobatto_IIIC_4_Solver_AC(u0, Mass, Stiff, f, f_prime, Tend, tN, F, al)

xN = length(Mass);
cut = ones(tN, 1);
energy = ones(tN+1,1);


tau = Tend/tN;
get_energy = @(u) (u' * Stiff * u) + ones(1, xN) * Mass * F(u);

energy(1) = get_energy(u0);
u(:,1) = u0;

Lap_m = Mass \ Stiff;

I = diag(ones(size(u0)));


RHS = @(u) - Lap_m * u + f(u);
RHS_prime = @(u) - Lap_m + diag(f_prime(u));
% a11 = 1/6;a12=-1/3;a13=1/6;a21=1/6;a22=5/12;a23=-1/12;a31=1/6;a32=2/3;a33=1/6;
sqrt5 = sqrt(5);
a11 = 1/12;a12=-sqrt5/12;a13=sqrt5/12;a14=-1/12;
a21 = 1/12;a22=1/4;a23=(10-7*sqrt5)/60;a24=sqrt5/60;
a31 = 1/12;a32=(10+7*sqrt5)/60;a33=1/4;a34=-sqrt5/60;
a41 = 1/12;a42=5/12;a43=5/12;a44=1/12;

b1=1/12;b2=5/12;b3=5/12;b4=1/12;

t1 = tic;
for i= 1:tN
    %         u(:,i+1) = Expl_RK4_Solver_AC(u(:,1),Mass, Stiff, f, tau, 100, F, al);
    
    ques = @(uu) [ ...
        u(:,i) + tau*a11*RHS(uu(1:end/4)) + tau*a12*RHS(uu((end/4+1):(2*end/4))) + tau*a13*RHS(uu((2*end/4+1):(3*end/4))) + tau*a14*RHS(uu((3*end/4+1):end)) - uu(1:end/4); ...
        u(:,i) + tau*a21*RHS(uu(1:end/4)) + tau*a22*RHS(uu((end/4+1):(2*end/4))) + tau*a23*RHS(uu((2*end/4+1):(3*end/4))) + tau*a24*RHS(uu((3*end/4+1):end)) - uu((end/4+1):(2*end/4)); ...
        u(:,i) + tau*a31*RHS(uu(1:end/4)) + tau*a32*RHS(uu((end/4+1):(2*end/4))) + tau*a33*RHS(uu((2*end/4+1):(3*end/4))) + tau*a34*RHS(uu((3*end/4+1):end)) - uu((2*end/4+1):(3*end/4)); ...
        u(:,i) + tau*a41*RHS(uu(1:end/4)) + tau*a42*RHS(uu((end/4+1):(2*end/4))) + tau*a43*RHS(uu((2*end/4+1):(3*end/4))) + tau*a44*RHS(uu((3*end/4+1):end)) - uu((3*end/4+1):end)] ;
    ques_prime = @(uu) [ ...
        tau*a11*RHS_prime(uu(1:end/4)) - I, tau*a12*RHS_prime(uu((end/4+1):(2*end/4))), tau*a13*RHS_prime(uu((2*end/4+1):(3*end/4))), tau*a14*RHS_prime(uu((3*end/4+1):end));...
        tau*a21*RHS_prime(uu(1:end/4)), tau*a22*RHS_prime(uu((end/4+1):(2*end/4))) - I, tau*a23*RHS_prime(uu((2*end/4+1):(3*end/4))), tau*a24*RHS_prime(uu((3*end/4+1):end));...
        tau*a31*RHS_prime(uu(1:end/4)), tau*a32*RHS_prime(uu((end/4+1):(2*end/4))), tau*a33*RHS_prime(uu((2*end/4+1):(3*end/4))) - I, tau*a34*RHS_prime(uu((3*end/4+1):end)); ...
        tau*a41*RHS_prime(uu(1:end/4)), tau*a42*RHS_prime(uu((end/4+1):(2*end/4))), tau*a43*RHS_prime(uu((2*end/4+1):(3*end/4))), tau*a44*RHS_prime(uu((3*end/4+1):end)) - I];
    x0 = [u(:,i);u(:,i);u(:,i);u(:,i)];
    answer = Newton_iter(ques, ques_prime, x0, 1e-12);
    
    u(:,i+1) = u(:,i) + tau*b1*RHS(answer(1:end/4)) + tau*b2*RHS(answer((end/4+1):(2*end/4))) + tau*b3*RHS(answer((2*end/4+1):(3*end/4))) + tau*b4*RHS(answer((3*end/4+1):end));
%     u(:,i+1) = fsolve(ques, u(:,i), opt);
%     u(:,i+1) = r1;
    
    cut(i) = max(max(abs(min(max(u(:,i+1),-al),al) - u(:,i+1)))); 
    u(:,i+1) = min(max(u(:,i+1),-al),al);  
    energy(i+1) = get_energy(u(:,i+1));
end

time1 = toc(t1);

uend = u(:,end);
end