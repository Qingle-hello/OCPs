function [uend, time1, cut, energy] = Lobatto_IIIC_2_Solver_AC(u0, Mass, Stiff, f, f_prime, Tend, tN, F, al)

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
a11 = 1/2;a12 = -1/2;a21 = 1/2;a22 = 1/2;
b1 = 1/2; b2 = 1/2;

t1 = tic;
for i= 1:tN
    %         u(:,i+1) = Expl_RK4_Solver_AC(u(:,1),Mass, Stiff, f, tau, 100, F, al);
    
    ques = @(uu) [u(:,i) + tau*a11*RHS(uu(1:end/2)) + tau*a12*RHS(uu((end/2+1):end)) - uu(1:end/2); ...
        u(:,i) + tau*a21*RHS(uu(1:end/2)) + tau*a22*RHS(uu((end/2+1):end)) - uu((end/2+1):end)] ;
    ques_prime = @(uu) [tau*a11*RHS_prime(uu(1:end/2)) - I, tau*a12*RHS_prime(uu((end/2+1):end));
        tau*a21*RHS_prime(uu(1:end/2)), tau*a22*RHS_prime(uu((end/2+1):end)) - I];
    x0 = [u(:,i);u(:,i)];
    answer = Newton_iter(ques, ques_prime, x0, 1e-12);
    
    u(:,i+1) = u(:,i) + tau*b1*RHS(answer(1:end/2)) + tau*b2*RHS(answer((end/2+1):end));
%     u(:,i+1) = fsolve(ques, u(:,i), opt);
%     u(:,i+1) = r1;
    
    cut(i) = max(max(abs(min(max(u(:,i+1),-al),al) - u(:,i+1)))); 
    u(:,i+1) = min(max(u(:,i+1),-al),al);  
    energy(i+1) = get_energy(u(:,i+1));
end
time1 = toc(t1);


uend = u(:,end);
end