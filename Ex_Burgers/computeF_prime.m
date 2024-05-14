function dfdu = computeF_prime(u, dx)
    N = length(u);
    dfdu = zeros(N, N); % Initialize a NxN matrix
    
    % 使用向量化方法计算 dfdu
    dfdu(2:N-1, 3:N) = diag(-u(3:N) / (2*dx));
    dfdu(2:N-1, 1:N-2) = dfdu(2:N-1, 1:N-2) + diag(u(1:N-2) / (2*dx));
end
