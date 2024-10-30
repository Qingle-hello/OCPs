function [uend, time1, cut, energy] = I_Euler_Solver_AC(u0, Mass, Stiff, f, f_prime, Tend, tN, F, al, L_fast, U_fast)

xN = length(Mass);
cut = ones(tN, 1);
energy = ones(tN+1,1);


tau = Tend/tN;
get_energy = @(u) (u' * Stiff * u) + ones(1, xN) * Mass * F(u);

energy(1) = get_energy(u0);
u(:,1) = u0;

Lap_m = Mass \ Stiff;

L = tau*Lap_m;
I = diag(ones(size(u0)));   
I_Euler_op = (I+L);    
% [L_IE,U_IE] = lu(I_Euler_op);
A = inv(I_Euler_op);

t1 = tic;
for i=1:tN

    f1 = f(u(:,i));
    u(:,i+1) = A*(u(:,i) + tau * f1 );     

    cut(i) = max(max(abs(min(max(u(:,i+1),-al),al) - u(:,i+1)))); 
    u(:,i+1) = min(max(u(:,i+1),-al),al);  
    energy(i+1) = get_energy(u(:,i+1));    

end    

time1 = toc(t1);

uend = u(:,end);

    

end