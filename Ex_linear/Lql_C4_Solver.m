function [uend,time] = Lql_C4_Solver(u0, tau, t_start, t_count, Mass,Stiff, f, x)

[NoF,~] = size(u0);

ulast = u0;
unext = u0;

Lap_m = Mass \ Stiff;

L = tau*Lap_m;
I = speye(NoF);

Lb3C2_op = (0.37898*L^2 + 0.78909*L + I);    
% [L_Lb3C2,U_Lb3C2] = lu(Lb3C2_op);
A = inv(Lb3C2_op);

tic;
for i=1:t_count 
    t2 = t_start + i * tau;
    f2 = f(t2, x);
    
    r = (0.37422*L + I) * f2;
    b = (0.98744*I - 0.22082*L)* ulast + tau * r;


    unext = A * b + 0.01256 * ulast;
    
   
    
    ulast = unext;
end    
time = toc;
uend = unext;

    

end
