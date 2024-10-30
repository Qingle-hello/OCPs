parpool('local', 60, 'IdleTimeout', 720000);

% coarse_solver = @ I_Euler_Solver_B;
% coarse_solver = @ Lql_C2_Solver_B;
% coarse_solver = @ Lql_C3_Solver_B;
% coarse_solver = @ Lql_CR_Solver_B;
% coarse_solver = @ Lql_C4_Solver_B;

% fine_solver = @ Radau_IIA_3_Solver_BI;
% fine_solver = @ Lobatto_IIIC_2_Solver_BI;
% fine_solver = @ Lobatto_IIIC_3_Solver_BI;
% fine_solver = @ Lobatto_IIIC_4_Solver_BI;

% nu = 1;  
% nu = 0.02;
% function Ex_3_a(coarse_solver, fine_solver, Nt_fine, J,nu, Nx)

% J = 10; Nt_fine = 2000; % nu = 1
% J =100; Nt_fine = 2000; % nu = 1
% J = 10; Nt_fine = 10000; % nu = 0.02
% J = 100; Nt_fine = 10000; % nu = 0.02

% Nx = 3000; 

%% Ex_3_a
Nx = 1000;
Ex_3_a(@ I_Euler_Solver_B,@ Lobatto_IIIC_2_Solver_BI, 20, 2000, 1,Nx)
Ex_3_a(@ I_Euler_Solver_B,@ Lobatto_IIIC_3_Solver_BI, 20, 2000, 1, Nx)
Ex_3_a(@ I_Euler_Solver_B,@ Radau_IIA_3_Solver_BI, 20, 2000, 1, Nx)
Ex_3_a(@ I_Euler_Solver_B,@ Lobatto_IIIC_4_Solver_BI, 20, 2000, 1, Nx)

Ex_3_a(@ Lql_C2_Solver_B, @ Lobatto_IIIC_2_Solver_BI, 20, 2000, 1,Nx)
Ex_3_a(@ Lql_C3_Solver_B, @ Lobatto_IIIC_3_Solver_BI, 20, 2000, 1,Nx)
Ex_3_a(@ Lql_CR_Solver_B, @ Radau_IIA_3_Solver_BI, 20, 2000, 1,Nx)
Ex_3_a(@ Lql_C4_Solver_B, @ Lobatto_IIIC_4_Solver_BI, 20, 2000, 1,Nx)

Ex_3_a(@ I_Euler_Solver_B,@ Lobatto_IIIC_2_Solver_BI, 100, 2000, 1,Nx)
Ex_3_a(@ I_Euler_Solver_B,@ Lobatto_IIIC_3_Solver_BI, 100, 2000, 1,Nx)
Ex_3_a(@ I_Euler_Solver_B,@ Radau_IIA_3_Solver_BI,  100, 2000,1, Nx)
Ex_3_a(@ I_Euler_Solver_B,@ Lobatto_IIIC_4_Solver_BI, 100, 2000, 1,Nx)

Ex_3_a(@ Lql_C2_Solver_B, @ Lobatto_IIIC_2_Solver_BI, 100, 2000, 1,Nx)
Ex_3_a(@ Lql_C3_Solver_B, @ Lobatto_IIIC_3_Solver_BI, 100, 2000, 1,Nx)
Ex_3_a(@ Lql_CR_Solver_B, @ Radau_IIA_3_Solver_BI, 100, 2000, 1,Nx)
Ex_3_a(@ Lql_C4_Solver_B, @ Lobatto_IIIC_4_Solver_BI, 100, 2000, 1,Nx)

%% Ex_3_b
Ex_3_b(@ I_Euler_Solver_B,@ Lobatto_IIIC_2_Solver_BI, 10, 10000, 0.02,Nx)
Ex_3_b(@ I_Euler_Solver_B,@ Lobatto_IIIC_3_Solver_BI, 10, 10000, 0.02,Nx)
Ex_3_b(@ I_Euler_Solver_B,@ Radau_IIA_3_Solver_BI, 10, 10000, 0.02,Nx)
Ex_3_b(@ I_Euler_Solver_B,@ Lobatto_IIIC_4_Solver_BI, 10, 10000, 0.02,Nx)

Ex_3_b(@ Lql_C2_Solver_B, @ Lobatto_IIIC_2_Solver_BI, 10, 10000, 0.02,Nx)
Ex_3_b(@ Lql_C3_Solver_B, @ Lobatto_IIIC_3_Solver_BI, 10, 10000, 0.02,Nx)
Ex_3_b(@ Lql_CR_Solver_B, @ Radau_IIA_3_Solver_BI, 10, 10000, 0.02,Nx)
Ex_3_b(@ Lql_C4_Solver_B, @ Lobatto_IIIC_4_Solver_BI, 10, 10000, 0.02,Nx)

Ex_3_b(@ I_Euler_Solver_B,@ Lobatto_IIIC_2_Solver_BI, 100, 10000, 0.02,Nx)
Ex_3_b(@ I_Euler_Solver_B,@ Lobatto_IIIC_3_Solver_BI, 100, 10000, 0.02,Nx)
Ex_3_b(@ I_Euler_Solver_B,@ Radau_IIA_3_Solver_BI, 100, 10000, 0.02,2*Nx)
Ex_3_b(@ I_Euler_Solver_B,@ Lobatto_IIIC_4_Solver_BI, 100, 10000, 0.02,2.5*Nx)

Ex_3_b(@ Lql_C2_Solver_B, @ Lobatto_IIIC_2_Solver_BI, 100, 10000, 0.02,Nx)
Ex_3_b(@ Lql_C3_Solver_B, @ Lobatto_IIIC_3_Solver_BI, 100, 10000, 0.02,Nx)
Ex_3_b(@ Lql_CR_Solver_B, @ Radau_IIA_3_Solver_BI, 100, 10000, 0.02,2*Nx)
Ex_3_b(@ Lql_C4_Solver_B, @ Lobatto_IIIC_4_Solver_BI, 100, 10000, 0.02,2.5*Nx)

Nx = 1000;
% J20 CR C4
Ex_3_a(@ Lql_CR_Solver_B, @ Radau_IIA_3_Solver_BI, 20, 2000, 1,Nx)
% Ex_3_a(@ Lql_C4_Solver_B, @ Lobatto_IIIC_4_Solver_BI, 20, 2000, 1,Nx)

Ex_3_a(@ Lql_CR_Solver_B, @ Radau_IIA_3_Solver_BI, 100, 2000, 1,Nx)
% Ex_3_a(@ Lql_C4_Solver_B, @ Lobatto_IIIC_4_Solver_BI, 100, 2000, 1,Nx)


fprintf('该程序已经结束');
pool = gcp('nocreate');  % 获取当前并行池，如果不存在则不创建
if ~isempty(pool)
    delete(pool);  % 如果并行池存在，删除并行池
end


