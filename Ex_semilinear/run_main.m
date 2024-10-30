% parpool('local', 60, 'IdleTimeout', 720000);

% coarse_solver = @ I_Euler_Solver_AC;
% coarse_solver = @ Lql_C2_Solver_AC;
% coarse_solver = @ Lql_C2_Solver_AC;
% coarse_solver = @ Lql_CR_Solver_AC;
% coarse_solver = @ Lql_C4_Solver_AC;

% fine_solver = @ Lobatto_IIIC_2_Solver_AC;
% fine_solver = @ Lobatto_IIIC_3_Solver_AC;
% fine_solver = @ Radau_IIA_3_Solver_AC;
% fine_solver = @ Lobatto_IIIC_4_Solver_AC;

% TN = 100; J = 20;  % fix tn = 2000, for ep = 1
% TN = 20; J = 100; % fix tn = 2000, for ep = 1
% TN = 1000; J = 10; % fix tn = 10000, for ep = sqrt(0.02)
% TN = 500; J = 20; % fix tn = 10000, for ep = sqrt(0.02)

% function Ex_2_b(coarse_solver, fine_solver, TN, J, ep)
%% Ex_2_a

% Ex_2_a(@ I_Euler_Solver_AC,@ Lobatto_IIIC_2_Solver_AC, 100, 20, 1)
% Ex_2_a(@ I_Euler_Solver_AC,@ Lobatto_IIIC_3_Solver_AC, 100, 20, 1)
% Ex_2_a(@ I_Euler_Solver_AC,@ Radau_IIA_3_Solver_AC, 100, 20, 1)
% Ex_2_a(@ I_Euler_Solver_AC,@ Lobatto_IIIC_4_Solver_AC, 100, 20, 1)
% 
% Ex_2_a(@ Lql_C2_Solver_AC, @ Lobatto_IIIC_2_Solver_AC, 100, 20, 1)
% Ex_2_a(@ Lql_C3_Solver_AC, @ Lobatto_IIIC_3_Solver_AC, 100, 20, 1)
% Ex_2_a(@ Lql_C4_Solver_AC, @ Radau_IIA_3_Solver_AC, 100, 20, 1)
% Ex_2_a(@ Lql_CR_Solver_AC, @ Lobatto_IIIC_4_Solver_AC, 100, 20, 1)

% Ex_2_a(@ I_Euler_Solver_AC,@ Lobatto_IIIC_2_Solver_AC, 20, 100, 1)
% Ex_2_a(@ I_Euler_Solver_AC,@ Lobatto_IIIC_3_Solver_AC, 20, 100, 1)
% Ex_2_a(@ I_Euler_Solver_AC,@ Radau_IIA_3_Solver_AC, 20, 100, 1)
% Ex_2_a(@ I_Euler_Solver_AC,@ Lobatto_IIIC_4_Solver_AC, 20, 100, 1)
% 
% Ex_2_a(@ Lql_C2_Solver_AC, @ Lobatto_IIIC_2_Solver_AC, 20, 100, 1)
% Ex_2_a(@ Lql_C3_Solver_AC, @ Lobatto_IIIC_3_Solver_AC, 20, 100, 1)
% Ex_2_a(@ Lql_C4_Solver_AC, @ Radau_IIA_3_Solver_AC, 20, 100, 1)
% Ex_2_a(@ Lql_CR_Solver_AC, @ Lobatto_IIIC_4_Solver_AC, 20, 100, 1)
% 
%% Ex_2_b
% Ex_2_b(@ I_Euler_Solver_AC,@ Lobatto_IIIC_2_Solver_AC, 1000, 10, 0.02)
% Ex_2_b(@ I_Euler_Solver_AC,@ Lobatto_IIIC_3_Solver_AC, 1000, 10, 0.02)
% Ex_2_b(@ I_Euler_Solver_AC,@ Radau_IIA_3_Solver_AC, 1000, 10, 0.02)
% Ex_2_b(@ I_Euler_Solver_AC,@ Lobatto_IIIC_4_Solver_AC, 1000, 10, 0.02)

% Ex_2_b(@ Lql_C2_Solver_AC, @ Lobatto_IIIC_2_Solver_AC, 1000, 10, 0.02)
% Ex_2_b(@ Lql_C3_Solver_AC, @ Lobatto_IIIC_3_Solver_AC, 1000, 10, 0.02)
% Ex_2_b(@ Lql_CR_Solver_AC, @ Radau_IIA_3_Solver_AC, 1000, 10, 0.02)
% Ex_2_b(@ Lql_CR_Solver_AC, @ Radau_IIA_3_Solver_AC, 1000, 10, 0.02)

% Ex_2_b(@ I_Euler_Solver_AC,@ Lobatto_IIIC_2_Solver_AC, 500, 20, 0.02)
% Ex_2_b(@ I_Euler_Solver_AC,@ Lobatto_IIIC_3_Solver_AC, 500, 20, 0.02)
% Ex_2_b(@ I_Euler_Solver_AC,@ Radau_IIA_3_Solver_AC, 500, 20, 0.02)
% Ex_2_b(@ I_Euler_Solver_AC,@ Lobatto_IIIC_4_Solver_AC, 500, 20, 0.02)

% Ex_2_b(@ Lql_C2_Solver_AC, @ Lobatto_IIIC_2_Solver_AC, 500, 20, 0.02)
% Ex_2_b(@ Lql_C3_Solver_AC, @ Lobatto_IIIC_3_Solver_AC, 500, 20, 0.02)
Ex_2_b(@ Lql_CR_Solver_AC, @ Radau_IIA_3_Solver_AC, 500, 20, 0.02)
% Ex_2_b(@ Lql_C4_Solver_AC, @ Lobatto_IIIC_4_Solver_AC, 500, 20, 0.02)


Ex_2_b(@ Lql_CR_Solver_AC, @ Radau_IIA_3_Solver_AC, 1000, 10, 0.02)
Ex_2_b(@ Lql_CR_Solver_AC, @ Radau_IIA_3_Solver_AC, 500, 20, 0.02)
fprintf('该程序已经结束');
pool = gcp('nocreate');  % 获取当前并行池，如果不存在则不创建
if ~isempty(pool)
    delete(pool);  % 如果并行池存在，删除并行池
end


