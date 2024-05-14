function [uend,time] = Radau_IIA_3_Solver(u0, tau, t_start, t_count, Mass, Stiff, f, x)

[N,~] = size(u0);

ulast = u0;
unext = u0;

Lap_m = Mass \ Stiff;

A = Lap_m;
I = speye(N);

L = -A;

s=sqrt(6);

g11 = (88-7*s)/360;g12=(296-169*s)/1800;g13=(-2+3*s)/225;
g21=(296+169*s)/1800;g22=(88+7*s)/360;g23=(-2-3*s)/225;
g31=(16-s)/36;g32=(16+s)/36;g33=1/9;

b1=(16-s)/36;b2=(16+s)/36;b3=1/9;


R = zeros(3*N);
R(1:N,1:N) = I + tau*g11*A;
R(1:N,N+1:2*N) = tau*g12*A;
R(1:N,2*N+1:3*N) = tau*g13*A;

R(N+1:2*N,1:N) = tau*g21*A;
R(N+1:2*N,N+1:2*N) = I + tau*g22*A;
R(N+1:2*N,2*N+1:3*N) = tau*g23*A;

R(2*N+1:3*N,1:N) = tau*g31*A;
R(2*N+1:3*N,N+1:2*N) = tau*g32*A;
R(2*N+1:3*N,2*N+1:3*N) = I + tau*g33*A;

G = zeros(3*N);
G(1:N,1:N) = g11*I;
G(1:N,N+1:2*N) = g12*I;
G(1:N,2*N+1:3*N) = g13*I;

G(N+1:2*N,1:N) = g21*I;
G(N+1:2*N,N+1:2*N) = g22*I;
G(N+1:2*N,2*N+1:3*N) = g23*I;

G(2*N+1:3*N,1:N) = g31*I;
G(2*N+1:3*N,N+1:2*N) = g32*I;
G(2*N+1:3*N,2*N+1:3*N) = g33*I;


F = zeros(3*N,1);
Un = zeros(3*N,1);
Un_last = zeros(3*N,1);

[L_R,U_R] = lu(R);

tic;
for i=1:t_count
    F(1:N,1) = f(t_start + (i+((-6-sqrt(6))/10)) * tau, x);
    F(N+1:2*N,1) = f(t_start + (i+((sqrt(6)-6)/10)) * tau, x);
    F(2*N+1:3*N,1) = f(t_start + i * tau, x);
    
    Un_last(1:N,1) = ulast;
    Un_last(N+1:2*N,1) = ulast;
    Un_last(2*N+1:3*N,1) = ulast;
    
    Un = (U_R\(L_R\(Un_last + tau*G*F)));
    
    unext = ulast + tau * (b1 * (L*Un(1:N,1)+F(1:N,1)) ...
            + b2 * (L*Un(N+1:2*N,1)+F(N+1:2*N,1))...
            + b3 * (L*Un(2*N+1:3*N,1)+F(2*N+1:3*N,1)));
           
    
    ulast = unext;
end    
time = toc;
uend = unext;

end
