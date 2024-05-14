function [uend,time] = Lql_C2_Solver(u0, tau, t_start, t_count, Mass,Stiff, f, x)

[NoF,~] = size(u0);

ulast = u0;
unext = u0;

Lap_m = Mass \ Stiff;

L = tau*Lap_m;
I = speye(NoF);

D = L^2;
Lb3C2_op = (0.37931*D + 0.79033*L + I);    
[L_Lb3C2,U_Lb3C2] = lu(Lb3C2_op);

tic;
for i=1:t_count 
    % t1 = t_start + (i-1) * tau;
    t2 = t_start + i * tau;
    f2 = f(t2, x);
    
    r = (0.37447*L + I) * f2;
    b = (I-0.20967*L + 0.00484*D)* ulast + tau * r;

    unext = (U_Lb3C2\(L_Lb3C2\(b)));    
    
    ulast = unext;
end    
time = toc;

uend = unext; 


end
