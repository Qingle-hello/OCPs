parpool('local', 3, 'IdleTimeout', 7200);

% coarse_solver = @ I_Euler_Solver;
% coarse_solver = @ Lql_C2_Solver;
% coarse_solver = @ Lql_C2_Solver;
% coarse_solver = @ Lql_C3_Solver;
% coarse_solver = @ Lql_C4_Solver;
% fine_solver = @ Radau_IIA_3_Solver;
% fine_solver = @ Lobatto_IIIC_2_Solver;
% fine_solver = @ Lobatto_IIIC_3_Solver;
% fine_solver = @ Lobatto_IIIC_4_Solver;

%% Ex_1_a
% Ex_1_a(@ I_Euler_Solver,@ Lobatto_IIIC_2_Solver, 20);
% Ex_1_a(@ I_Euler_Solver,@ Lobatto_IIIC_3_Solver, 20);
% Ex_1_a(@ I_Euler_Solver,@ Lobatto_IIIC_4_Solver, 20);
% Ex_1_a(@ I_Euler_Solver,@ Radau_IIA_3_Solver, 20);
% 
% Ex_1_a(@ Lql_C2_Solver,@ Lobatto_IIIC_2_Solver, 20);
% Ex_1_a(@ Lql_C3_Solver,@ Lobatto_IIIC_3_Solver, 20);
% Ex_1_a(@ Lql_C4_Solver,@ Lobatto_IIIC_4_Solver, 20);
% Ex_1_a(@ Lql_CR_Solver,@ Radau_IIA_3_Solver, 20);

Ex_1_a(@ I_Euler_Solver,@ Lobatto_IIIC_2_Solver, 100);
Ex_1_a(@ I_Euler_Solver,@ Lobatto_IIIC_3_Solver, 100);
Ex_1_a(@ I_Euler_Solver,@ Lobatto_IIIC_4_Solver, 100);
Ex_1_a(@ I_Euler_Solver,@ Radau_IIA_3_Solver, 100);

Ex_1_a(@ Lql_C2_Solver,@ Lobatto_IIIC_2_Solver, 100);
Ex_1_a(@ Lql_C3_Solver,@ Lobatto_IIIC_3_Solver, 100);
Ex_1_a(@ Lql_C4_Solver,@ Lobatto_IIIC_4_Solver, 100);
Ex_1_a(@ Lql_CR_Solver,@ Radau_IIA_3_Solver, 100);

%% Ex_1_b
% Ex_1_b(@ I_Euler_Solver,@ Lobatto_IIIC_2_Solver, 20);
% Ex_1_b(@ I_Euler_Solver,@ Lobatto_IIIC_3_Solver, 20);
% Ex_1_b(@ I_Euler_Solver,@ Lobatto_IIIC_4_Solver, 20);
% Ex_1_b(@ I_Euler_Solver,@ Radau_IIA_3_Solver, 20);
% 
% Ex_1_b(@ Lql_C2_Solver,@ Lobatto_IIIC_2_Solver, 20);
% Ex_1_b(@ Lql_C3_Solver,@ Lobatto_IIIC_3_Solver, 20);
% Ex_1_b(@ Lql_C4_Solver,@ Lobatto_IIIC_4_Solver, 20);
% Ex_1_b(@ Lql_CR_Solver,@ Radau_IIA_3_Solver, 20);
% 
% Ex_1_b(@ I_Euler_Solver,@ Lobatto_IIIC_2_Solver, 100);
% Ex_1_b(@ I_Euler_Solver,@ Lobatto_IIIC_3_Solver, 100);
% Ex_1_b(@ I_Euler_Solver,@ Lobatto_IIIC_4_Solver, 100);
% Ex_1_b(@ I_Euler_Solver,@ Radau_IIA_3_Solver, 100);
% 
% Ex_1_b(@ Lql_C2_Solver,@ Lobatto_IIIC_2_Solver, 100);
% Ex_1_b(@ Lql_C3_Solver,@ Lobatto_IIIC_3_Solver, 100);
% Ex_1_b(@ Lql_C4_Solver,@ Lobatto_IIIC_4_Solver, 100);
% Ex_1_b(@ Lql_CR_Solver,@ Radau_IIA_3_Solver, 100);


fprintf('该程序已经结束');
pool = gcp('nocreate');  % 获取当前并行池，如果不存在则不创建
if ~isempty(pool)
    delete(pool);  % 如果并行池存在，删除并行池
end


