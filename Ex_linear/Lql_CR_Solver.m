function [uend,time] = Lql_CR_Solver(u0, tau, t_start, t_count, Mass, Stiff, f, x)

[NoF,~] = size(u0);

ulast = u0;
unext = u0;

Lap_m = Mass \ Stiff;

L = tau*Lap_m;
I = speye(NoF);
Lb3C2_op = (0.37922*L^2 + 0.78875*L + I);    
[L_Lb3C2,U_Lb3C2] = lu(Lb3C2_op);

tic;
for i=1:t_count 
    t2 = t_start + i * tau;
    f2 = f(t2, x);
    
    r = (0.37443*L + I) * f2;
    b = (I-0.21125*L + 0.00479*L^2)* ulast + tau * r;


    unext = (U_Lb3C2\(L_Lb3C2\(b)));    
    
   
    
    ulast = unext;
end    
time = toc;
uend = unext;

    

end
