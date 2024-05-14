function [uend,time] = Lobatto_IIIC_4_Solver(u0, tau, t_start, t_count, Mass, Stiff, f, x)

[N,~] = size(u0);

ulast = u0;
unext = u0;

Lap_m = Mass \ Stiff;
A = Lap_m; % positive defined
I = speye(N);

L = -A;


sqrt5 = sqrt(5);
g11 = 1/12;g12=-sqrt5/12;g13=sqrt5/12;g14=-1/12;
g21 = 1/12;g22=1/4;g23=(10-7*sqrt5)/60;g24=sqrt5/60;
g31 = 1/12;g32=(10+7*sqrt5)/60;g33=1/4;g34=-sqrt5/60;
g41 = 1/12;g42=5/12;g43=5/12;g44=1/12;

R = zeros(4*N);
R(1:N,1:N) = I + tau*g11*A;
R(1:N,N+1:2*N) = tau*g12*A;
R(1:N,2*N+1:3*N) = tau*g13*A;
R(1:N,3*N+1:4*N) = tau*g14*A;

R(N+1:2*N,1:N) = tau*g21*A;
R(N+1:2*N,N+1:2*N) = I + tau*g22*A;
R(N+1:2*N,2*N+1:3*N) = tau*g23*A;
R(N+1:2*N,3*N+1:4*N) = tau*g24*A;

R(2*N+1:3*N,1:N) = tau*g31*A;
R(2*N+1:3*N,N+1:2*N) = tau*g32*A;
R(2*N+1:3*N,2*N+1:3*N) = I + tau*g33*A;
R(2*N+1:3*N,3*N+1:4*N) = tau*g34*A;

R(3*N+1:4*N,1:N) = tau*g41*A;
R(3*N+1:4*N,N+1:2*N) = tau*g42*A;
R(3*N+1:4*N,2*N+1:3*N) = tau*g43*A;
R(3*N+1:4*N,3*N+1:4*N) = I + tau*g44*A;

G = zeros(4*N);
G(1:N,1:N) = g11*I;
G(1:N,N+1:2*N) = g12*I;
G(1:N,2*N+1:3*N) = g13*I;
G(1:N,3*N+1:4*N) = g14*I;

G(N+1:2*N,1:N) = g21*I;
G(N+1:2*N,N+1:2*N) = g22*I;
G(N+1:2*N,2*N+1:3*N) = g23*I;
G(N+1:2*N,3*N+1:4*N) = g24*I;

G(2*N+1:3*N,1:N) = g31*I;
G(2*N+1:3*N,N+1:2*N) = g32*I;
G(2*N+1:3*N,2*N+1:3*N) = g33*I;
G(2*N+1:3*N,3*N+1:4*N) = g34*I;

G(3*N+1:4*N,1:N) = g41*I;
G(3*N+1:4*N,N+1:2*N) = g42*I;
G(3*N+1:4*N,2*N+1:3*N) = g43*I;
G(3*N+1:4*N,3*N+1:4*N) = g44*I;

F = zeros(4*N,1);
Un = zeros(4*N,1);
Un_last = zeros(4*N,1);

[L_R,U_R] = lu(R);

tic;

for i=1:t_count
    F(1:N,1) = f(t_start + (i-1) * tau, x);
    F(N+1:2*N,1) = f(t_start + (i-0.5-sqrt(5)/10) * tau, x);
    F(2*N+1:3*N,1) = f(t_start + (i-0.5+sqrt(5)/10) * tau, x);
    F(3*N+1:4*N,1) = f(t_start + i * tau, x);
    
    Un_last(1:N,1) = ulast;
    Un_last(N+1:2*N,1) = ulast;
    Un_last(2*N+1:3*N,1) = ulast;
    Un_last(3*N+1:4*N,1) = ulast;
    
    Un = (U_R\(L_R\(Un_last + tau*G*F)));
    
    unext = ulast + tau * (1/12 * (L*Un(1:N,1)+F(1:N,1)) ...
            + 5/12 * (L*Un(N+1:2*N,1)+F(N+1:2*N,1))...
            + 5/12 * (L*Un(2*N+1:3*N,1)+F(2*N+1:3*N,1))...
            + 1/12 * (L*Un(3*N+1:4*N,1)+F(3*N+1:4*N,1)));
    
    ulast = unext;
end    
time = toc;
uend = unext;

end
