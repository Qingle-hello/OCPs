function [uend,time] = Lql_C3_Solver(u0, tau, t_start, t_count, Mass,Stiff ,f, x)

[NoF,~] = size(u0);

ulast = u0;
unext = u0;

Lap_m = Mass \ Stiff;

L = tau*Lap_m;
I = speye(NoF);
D=L*L;
Lb3C3_op = I+0.78986*L+0.38283*D;   
[L_Lb3C3,U_Lb3C3] = lu(Lb3C3_op);

tic;
for i=1:t_count
    f2 = f(t_start + i * tau, x);
    r = (0.37797*L+I) * f2;
    r = (I-0.21014*L+0.00486*D)*ulast + tau * r;
    unext = (U_Lb3C3\(L_Lb3C3\r));    

    ulast = unext;
end    

time = toc;
uend = unext;

    

end
