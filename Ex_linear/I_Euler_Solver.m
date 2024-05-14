function [uend,time] = I_Euler_Solver(u0, tau, t_start, t_count, Mass,Stiff, f, x)

[NoF,~] = size(u0);

ulast = u0;
unext = u0;

Lap_m = Mass \ Stiff;

L = tau*Lap_m;
I = speye(NoF);
I_Euler_op = (I+L);    
[L_IE,U_IE] = lu(I_Euler_op);
tic;

for i=1:t_count
    fn = f(t_start + i * tau, x);
    unext = (U_IE\(L_IE\(ulast + tau * fn )));    
        
    ulast = unext;
end 
time = toc;   

uend = unext;

    

end
