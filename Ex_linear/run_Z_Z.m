clc;clear;close all;

L = pi;      % inteval x \in [0, L];
M0 = 100;   % count of sub-intervals of space
r = 1;         % order on space

T = 1;     % t \in [0, T];
tn = 100;  % count of sub-intervals of time
             
J = 10;       % mesh ratio
iter_max = 30;

% 记得改save
% Lobatto_IIIC_2_Solver
% Lobatto_IIIC_3_Solver
% Exact_Solver TODO
% coarse_solver = @ I_Euler_Solver;
coarse_solver = @ Lql_CR_Solver;
% coarse_solver = @ CN_Solver;
% fine_solver = @ Lobatto_IIIC_3_Solver;
% fine_solver = @ Lobatto_IIIC_3_Solver;
fine_solver = @ Radau_IIA_3_Solver;

% f = @(t, x) 50*sin(2*pi*(x+t));
f = @(t, x) cos(t) .* sin(x);

% get_u0 = @(x)  ((x.^10) .* ( x-pi).^10)./((pi/2).^10) ;
get_u0 = @(x) (  (x<=pi/2).*1 + (x>pi/2).*(-1)  );
% get_u0 = @(x) sin(x);
% u_t = Laplacian u + f(t) 
% get_exact_u = @(t, x) exp(-t) * cos(x);
% get_exact_u = @(t, x)  cos(x);

TN = tn / J;
TAU = T / TN;
tau = T / tn;

[Mass, Stiff, x0, NoF]=fem_matrix_nm(L,M0,r); %Stiff=Stiff;
% Stiff = Mass + Stiff; % u_t - Lap u + u = f;

u0 = get_u0(x0);
time = zeros(2,TN);
total_time_all = zeros(iter_max,TN);

[u_seq_end,time_fine] = fine_solver(u0, tau, 0, tn, Mass, Stiff, f, x0);  

u_exact_end = u_seq_end;
% save('parareal efficiency/J_100_Exact_Solver.mat','u_exact_end');
% load('parareal efficiency/J_100_Exact_Solver.mat');

% u_err_para_seq = u_exact_end - u_seq_end;
% err(1) = sqrt(u_err_para_seq' * Mass * u_err_para_seq);
% disp(err(1));

% save('parareal efficiency 2/3-stage Radau/exact_time.mat','time_fine');

u_para(:,1) = u0;
u_para_new(:,1) = u0;

total_time = 0; % Initialize a variable to hold the total elapsed time

for k = 1 : TN
    [u_para(:, k+1),time(1,k)] = coarse_solver(u_para(:, k), TAU, (k-1)*TAU, 1, Mass, Stiff, f, x0);
    % u_para(:, k+1) = u0;
end

u_err_para_seq = u_para(:,end) - u_exact_end;
err(1) = sqrt(u_err_para_seq' * Mass * u_err_para_seq);

u_temp = zeros(M0-1,TN+1);

num_workers = min(TN,6);
parpool(num_workers);

for m = 1:iter_max
    parfor k = 1 : TN
        [u_temp(:,k+1),time(2,k)] = fine_solver(u_para(:, k), tau, (k-1)*TAU, J, Mass, Stiff, f, x0);
    end
    total_time_all(m,:) = time(2,:);

    t1 = tic;    
    for k = 1 : TN
        u_para_new(:, k+1) = ...
            coarse_solver(u_para_new(:, k), TAU, (k-1)*TAU, 1, Mass, Stiff, f, x0) ...
            - coarse_solver(u_para(:, k), TAU, (k-1)*TAU, 1, Mass, Stiff, f, x0) ...
            + u_temp(:,k+1);
    end
    u_para = u_para_new;
    
    u_err_para_seq = u_para(:,end) - u_exact_end;
    err(m+1) = sqrt(u_err_para_seq' * Mass * u_err_para_seq);

    % Stop the timer for this iteration and add to total_time
    elapsed_time = toc(t1);
    total_time = total_time + elapsed_time;
    
    % Print out the elapsed time for this iteration and the total time so far
    fprintf('Elapsed time for iteration %d: %.2f seconds\n', m, elapsed_time);
    fprintf('Total elapsed time: %.2f seconds\n', total_time);
    
end

Exact_time = 269;
time_fine = sum(max(total_time_all, [], 2));
speed_up_1 = (Exact_time)/(m * sum(time(1,:)) + time_fine);

eff_1 = 100*speed_up_1/(TN);

speed_up_2 = (Exact_time)/(time_fine);
eff_2 = 100*speed_up_2/(TN);

% save('parareal efficiency 2/4-stage Lobatto/iteration.mat','m');
% save('parareal efficiency 2/4-stage Lobatto/speed up with cost.mat','speed_up_1');
% save('parareal efficiency 2/4-stage Lobatto/efficiency with cost.mat','eff_1');

% save('parareal efficiency 2/4-stage Lobatto/speed up without cost.mat','speed_up_2');
% save('parareal efficiency 2/4-stage Lobatto/efficiency without cost.mat','eff_2');

% save('parareal efficiency 2/4-stage Lobatto/time.mat','total_time_all');

% save('parareal efficiency 2/4 BE/iteration.mat','m');
% save('parareal efficiency 2/4 BE/speed up with cost.mat','speed_up_1');
% save('parareal efficiency 2/4 BE/efficiency with cost.mat','eff_1');

% save('parareal efficiency 2/4 BE/speed up without cost.mat','speed_up_2');
% save('parareal efficiency 2/4 BE/efficiency without cost.mat','eff_2');

% save('parareal efficiency 2/4 BE/time.mat','total_time_all');

delete(gcp('nocreate'));

m_list = 0:iter_max; 

% Plot Figure 1 (Semilogy plot)
figure(1);
semilogy(m_list, err, '-o', 'DisplayName', 'Error');
hold on;
ref1 = (0.3).^(m_list - 1);
semilogy(m_list, ref1, '--', 'DisplayName', '0.3^k');
hold off;
xlabel('Iterations');
ylabel('Error (log scale)');
legend('Location', 'best');
title('Coarse propagator:q=1, fine propagator Labatto_IIIC_3_Solver');

% save('smooth/J100/J_100_Radau_IIA_3_Solver.mat', 'err');
% save('smooth/J100/J_100_Lobatto_IIIC_4_Solver.mat', 'err');
% save('smooth/J100/J_100_Lql_C4_Solver.mat', 'err');

% save('smooth/J20/J_20_Radau_IIA_3_Solver.mat', 'err');
% save('smooth/J20/J_20_Lobatto_IIIC_4_Solver.mat', 'err');
% save('smooth/J20/J_20_Lql_C4_Solver.mat', 'err');


% Save the 'err' variable to a file in the 'plot' folder
% save('robust/J100/J_100_Radau_IIA_3_Solver.mat', 'err');
% save('robust/J100/J_100_Lobatto_IIIC_4_Solver.mat', 'err');
% save('robust/J100/J_100_Lql_C4_Solver.mat', 'err');

% save('robust/J20/J_20_Radau_IIA_3_Solver.mat', 'err');
% save('robust/J20/J_20_Lobatto_IIIC_4_Solver.mat', 'err');
% save('robust/J20/J_20_Lql_C4_Solver.mat', 'err');
%% 数值收敛阶
alpha1=zeros(iter_max-1, 1);
for i=1:1:iter_max-1
    alpha1(i)=exp((log(err(i+1))-log(err(i))));
end

