function [x] = Newton_iter(f, f_prime, x0, tol)
if nargin < 4
    tol = 1e-12;
end
x_cur = x0;
lasterr = 2 * norm(f(x_cur));
thiserr = norm(f(x_cur));
maxtrial = 50;
trial = 0;
while (thiserr >= tol) && (thiserr <= 0.5 * lasterr) && (trial < maxtrial)
    lasterr = thiserr;
    fp = f_prime(x_cur);
    x_cur = x_cur - fp \ f(x_cur);
    thiserr = norm(f(x_cur));
    trial = trial + 1;
end
x = x_cur;
if (trial >= maxtrial)
    disp("尝试次数过多")
end
if (thiserr > 0.5 * lasterr) && (thiserr > 1e-9)
    disp("可能不收敛")
end
end
