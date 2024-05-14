% Example 4.1, problem data (a):
% T = 1,u0 = x^10(x-pi)^10 / (pi/2)^10, f = 0;

clc;clear;close all;

L = pi;      % inteval x \in [0, L];
M0 = 1000;   % count of sub-intervals of space
r = 1;         % order on space

T = 1;     % t \in [0, T];

tn = 2000;  % count of sub-intervals of fine time step
% J = 20;       % mesh ratio
J = 100;

iter_max = 10;

% coarse_solver = @ I_Euler_Solver;
coarse_solver = @ Lql_C2_Solver;
% coarse_solver = @ Lql_C2_Solver;
% coarse_solver = @ Lql_C3_Solver;
% coarse_solver = @ Lql_C4_Solver;

% fine_solver = @ Radau_IIA_3_Solver;
fine_solver = @ Lobatto_IIIC_2_Solver;
% fine_solver = @ Lobatto_IIIC_3_Solver;
% fine_solver = @ Lobatto_IIIC_4_Solver;

% f = @(t, x) 50*sin(2*pi*(x+t));
% f = @(t, x) cos(t) .* sin(x);
f = @(t, x) 0.*x;

get_u0 = @(x)  ((x.^10) .* ( x-pi).^10)./((pi/2).^10) ;
% get_u0 = @(x) (  (x<=pi/2).*1 + (x>pi/2).*(-1)  );


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


for m = 1:iter_max
    t1 = tic; 
    % parfor
    for k = 1 : TN
        [u_temp(:,k+1),time(2,k)] = fine_solver(u_para(:, k), tau, (k-1)*TAU, J, Mass, Stiff, f, x0);
    end
    total_time_all(m,:) = time(2,:);

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

%% convergence rate
alpha1=zeros(iter_max-1, 1);
for i=1:1:iter_max-1
    alpha1(i)=exp((log(err(i+1))-log(err(i))));
end

fprintf('The convergence rate: \n');
disp(alpha1);


m_list = 0:iter_max; 

% Plot Figure 1 (Semilogy plot)
figure(1);
semilogy(m_list, err, '-o', 'DisplayName', 'Error');
hold on;
% ref1 = (0.3).^(m_list - 1);
ref1 = (0.014).^(m_list - 1);
semilogy(m_list, ref1, '--', 'DisplayName', '0.014^k');
hold off;
xlabel('Iterations');
ylabel('Error (log scale)');
legend('Location', 'best');
title('LCP');


