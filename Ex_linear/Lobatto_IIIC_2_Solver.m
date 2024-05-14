function [uend,time] = Lobatto_IIIC_2_Solver(u0, tau, t_start, t_count, Mass,Stiff, f, x)

[N,~] = size(u0);

ulast = u0;
unext = u0;

Lap_m = Mass \ Stiff;

A = Lap_m;
I = speye(N);

L = -A;

g11 = 1/2;g12 = -1/2;g21 = 1/2;g22 = 1/2;
b1 = 1/2; b2 = 1/2;

R = zeros(2*N);
R(1:N,1:N) = I + tau*g11*A;
R(1:N,N+1:2*N) = tau*g12*A;

R(N+1:2*N,1:N) = tau*g21*A;
R(N+1:2*N,N+1:2*N) = I + tau*g22*A;


G = zeros(2*N);
G(1:N,1:N) = g11*I;
G(1:N,N+1:2*N) = g12*I;

G(N+1:2*N,1:N) = g21*I;
G(N+1:2*N,N+1:2*N) = g22*I;


F = zeros(2*N,1);
Un = zeros(2*N,1);
Un_last = zeros(2*N,1);

[L_R,U_R] = lu(R);

time1 = tic;
for i=1:t_count
    t1 = t_start + (i-1) * tau; 
    t2 = t_start + i * tau;
    F(1:N,1) = f(t1, x);
    F(N+1:2*N,1) = f(t2, x);

    Un_last(1:N,1) = ulast;
    Un_last(N+1:2*N,1) = ulast;
    
    Un = (U_R\(L_R\(Un_last + tau*G*F)));
    
    unext = ulast + tau * (b1 * (L*Un(1:N,1)+F(1:N,1)) ...
            + b2 * (L*Un(N+1:2*N,1)+F(N+1:2*N,1)));
    
    ulast = unext;
end    
time = toc(time1);
uend = unext;

end
