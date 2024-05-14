function fUpwind = computeF(u, dx)
    N = length(u);
    fUpwind = zeros(N, 1);
    
    % 使用向量化方法计算 fUpwind
    fUpwind(2:N-1) = -(u(3:N).^2 - u(1:N-2).^2) / (4*dx);
end
