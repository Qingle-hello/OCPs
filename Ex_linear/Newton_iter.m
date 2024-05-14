function [x] = Newton_iter(f, f_prime, x0, tol)
if nargin < 4
    tol = 1e-12;
end
x_cur = x0;
while norm(f(x_cur)) >= tol
    fp = f_prime(x_cur);
    x_cur = x_cur - fp \ f(x_cur);
end
x = x_cur;
end
