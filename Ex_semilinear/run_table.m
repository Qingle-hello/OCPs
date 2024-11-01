clc;clear;close all;

tic;
L=pi;       % inteval x \in [0, L];
M0 = 1000;    % count of sub-intervals of space
r = 1;        % order on space
T = 10;      % t \in [0, T];

% TN is the number of coarse steps
% here the target accuracy is 1e-6
TN = 150; J = 20; % three-stage Lobatto
% TN = 50; J = 20; % three-stage Radau
% TN = 40; J = 20; % four-stage Lobatto

iter_max = min([TN, 10]);

ep = sqrt(1);

% coarse_solver = @ I_Euler_Solver_AC;
coarse_solver = @ Lql_C2_Solver_AC;
% coarse_solver = @ Lql_C3_Solver_AC;
% coarse_solver = @ Lql_CR_Solver_AC;
% coarse_solver = @ Lql_C4_Solver_AC;

fine_solver = @ Lobatto_IIIC_2_Solver_AC;
% fine_solver = @ Lobatto_IIIC_3_Solver_AC;
% fine_solver = @ Radau_IIA_3_Solver_AC;
% fine_solver = @ Lobatto_IIIC_4_Solver_AC;


f = @(u) (u-u.^3)/(ep*ep);
f_prime = @(u) ((u./u) - 3 * u.^2) / (ep*ep);
F = @(u) (u.^2-1).^2/(4*ep*ep);
% f = @(u) 0*u;
% F = @(u) 0*u;

al = 1;

get_u0 = @(x) (  (x<=pi/2).*1 + (x>pi/2).*(-1));


tn = TN * J;
TAU = T / TN;
tau = T / tn;

[Mass, Stiff, x0, NoF]=fem_matrix_nm(L,M0,r); %Stiff=Stiff;
% Stiff = Mass + Stiff; % u_t - Lap u + u = f;

u0 = get_u0(x0);

%%fast IE on coarse
Lap_m = Mass \ Stiff;

L = TAU*Lap_m;
I = diag(ones(size(u0)));   
I_Euler_op = (I+L);    
[L_c,U_c] = lu(I_Euler_op);

time_total = zeros(iter_max+1,TN);

u_seq(:,1) = u0;
u_para(:,1) = u0;
u_para_new(:,1) = u0;

for k = 1 : TN
    fprintf("sequential solution: k=%d\t%f\n", k, toc)
    u_seq(:,k+1) = fine_solver(u_seq(:,k), Mass, Stiff, f, f_prime, TAU, J, F, al);
    % initialization
    [u_para(:, k+1),time_total(1,k)] = coarse_solver(u_para(:, k), Mass, Stiff, f, f_prime, TAU, 1, F, al, L_c,U_c);
    % u_para(:, k+1) = u0;
end

load('parareal efficiency/J_20_Exact_Solver.mat');

u_err_para_seq = u_para(:,2:end) - u_exact_end;
this_err = zeros(size(1:TN));
for k2 = 1 : TN
    this_err(k2) = sqrt(u_err_para_seq(:,k2)' * Mass * u_err_para_seq(:,k2));
end

err(1) = max(this_err);

%% parareal

time_temp = zeros(1,TN);
for m = 1:iter_max
    iteration_start_time = tic;
    % parfor k = 1 : TN
    for k = 1 : TN
        u_temp(:,k+1) = fine_solver(u_para(:, k), Mass, Stiff, f, f_prime, TAU, J, F, al);
    end

    for k = 1 : TN
        u_para_new(:, k+1) = ...
            coarse_solver(u_para_new(:, k), Mass, Stiff, f, f_prime, TAU, 1, F, al, L_c,U_c) ...
            - coarse_solver(u_para(:, k), Mass, Stiff, f, f_prime, TAU, 1, F, al, L_c,U_c) ...
            + u_temp(:,k+1);
    end
    u_para = u_para_new;
    
    u_err_para_seq = u_para(:,2:end) - u_exact_end;
    this_err = zeros(size(1:TN));
    for k2 = 1 : TN
        this_err(k2) = sqrt(u_err_para_seq(:,k2)' * Mass * u_err_para_seq(:,k2));
    end
    err(m+1) = max(this_err);

    fprintf("Parallel iteration %d took %f seconds, %f minutes more \n", m, toc(iteration_start_time),toc(iteration_start_time)*(iter_max-m)/60);
end

 
% Save the 'err' variable to a file in the 'plot' folder
% save('new J10/J_10_Radau_IIA_3_Solver.mat', 'err');
% save('new J10/J_10_Lobatto_IIIC_2_Solver.mat', 'err');
% save('new J10/J_10_Lql_C2_Solver.mat', 'err');

% save('new J20/J_20_Radau_IIA_3_Solver.mat', 'err');
% save('new J20/J_20_Lobatto_IIIC_2_Solver.mat', 'err');
% save('new J20/J_20_Lql_C2_Solver.mat', 'err');


%% convergence rate
alpha1=zeros(iter_max-1, 1);
for i=1:1:iter_max-1
    alpha1(i)=exp((log(err(i+1))-log(err(i))));
end

disp(alpha1);

%% Parareal efficiency
Exact_time = 2168; % three-stage Lobatto
% Exact_time = 597; % three-stage Radau
% Exact_time = 888; % four-stage Lobatto

% time_fine = sum(max(time_total(2:end,:), [], 2));

% speed_up_1 = (Exact_time)/(m * sum(time_total(1,:)) + time_fine);
% 
% eff_1 = 100*speed_up_1/(TN);
% 
% speed_up_2 = (Exact_time)/(time_fine);
% eff_2 = 100*speed_up_2/(TN);

% save('parareal efficiency 2/2-stage Lobatto/iteration.mat','m');
% save('parareal efficiency 2/2-stage Lobatto/speed up with cost.mat','speed_up_1');
% save('parareal efficiency 2/2-stage Lobatto/efficiency with cost.mat','eff_1');
% save('parareal efficiency 2/2-stage Lobatto/speed up without cost.mat','speed_up_2');
% save('parareal efficiency 2/2-stage Lobatto/efficiency without cost.mat','eff_2');

% save('parareal efficiency 2/2 BE/iteration.mat','m');
% save('parareal efficiency 2/2 BE/speed up with cost.mat','speed_up_1');
% save('parareal efficiency 2/2 BE/efficiency with cost.mat','eff_1');
% 
% save('parareal efficiency 2/2 BE/speed up without cost.mat','speed_up_2');
% save('parareal efficiency 2/2 BE/efficiency without cost.mat','eff_2');

