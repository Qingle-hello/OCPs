clc;clear;close all;

% Burgers' equation parameters and initial setup
Nx = 3000; 
L = 1;
T = 1;

nu = 1;  
% nu = 0.02;

J = 10; Nt_fine = 2000; % nu = 1
% J =100; Nt_fine = 2000; % nu = 1
% J = 10; Nt_fine = 10000; % nu = 0.02
% J = 100; Nt_fine = 10000; % nu = 0.02

dx = L / Nx;
iter_max = 10;
% Define coarse and fine time steps
Nt_coarse = Nt_fine/J;

dt_coarse = T/Nt_coarse; 
dt_fine = T/Nt_fine;


x = linspace(0, L, Nx+1)';
% u0 = sin(2*pi*x) ;

% u0 = sin(2*pi*x) + 0.5*sin(pi*x);
u0 =   (x<=1/2).*1 + (x>1/2).*(0) ;

u0(1) = 0; u0(end) = 0; % Boundary conditions
u = u0;
unew = u0;


% coarse_solver = @ I_Euler_Solver_B;
coarse_solver = @ Lql_C2_Solver_B;
% coarse_solver = @ Lql_C3_Solver_B;
% coarse_solver = @ Lql_CR_Solver_B;
% coarse_solver = @ Lql_C4_Solver_B;

% fine_solver = @ Radau_IIA_3_Solver_BI;
fine_solver = @ Lobatto_IIIC_2_Solver_BI;
% fine_solver = @ Lobatto_IIIC_3_Solver_BI;
% fine_solver = @ Lobatto_IIIC_4_Solver_BI;

% Matrix to store u values for all times
u_para = zeros(Nx+1, Nt_coarse+1);
u_seq = zeros(Nx+1,Nt_coarse+1);  % exact solution
u_temp = zeros(Nx+1,Nt_coarse+1); % record fine 
err = zeros(iter_max+1,1);


% 创建一个图形窗口
h = figure;

A = computeA(Nx,nu,dx);

u_seq(:,1) = u0;
u_para(:,1) = u0;
u_para_new(:,1) = u0;
u_temp(:,1) = u0;

%% coarse initialization
for k = 1 : Nt_coarse
    iteration_start_time = tic;
    u_seq(:,k+1) = fine_solver(u_seq(:,k),A,dx,dt_fine,J);
    u_para(:, k+1) = coarse_solver(u_para(:,k),A,dx,dt_coarse);
    % u_para(:, k+1) = u0; % for plot
    fprintf("Iteration %d took %f seconds. %f minutes more\n", k, toc(iteration_start_time),toc(iteration_start_time)*(Nt_coarse-k)/60);
end

u_seq_end = u_seq(:,end);

u_err_para_seq = u_para(:,2:end) - u_seq(:,2:end);
    this_err = zeros(size(1:Nt_coarse));
    for k2 = 1 : Nt_coarse
        this_err(k2) = sqrt(1/Nx * u_err_para_seq(:,k2)' * u_err_para_seq(:,k2));
    end
    
err(1) = max(this_err);

% Flipping u_para by rows
    u_para_flipped = flipud(u_seq);

    % Plot the figure
    h_surf = surf(linspace(0, T, Nt_coarse+1), x, u_para_flipped, 'EdgeColor', 'none');
    xlabel('Time (t)');
    ylabel('Position (x)');
    zlabel('u(x,t)');
   % title(['u(x,t) over time and position, m = ' num2str(m) ', err = ' num2str(err(m))]);
    view(3); % 3D view

    % Set colormap and increase number of levels
    % colormap(jet(256)); % 'jet' colormap with 256 levels
    colorbar;


    Yticks_locs = get(gca, 'YTick');

    Yticks_labels = L - Yticks_locs;

    set(gca, 'YTickLabel', arrayfun(@num2str, Yticks_labels, 'UniformOutput', false));

    pause(0.05);

    

%% parareal
for m = 1:iter_max
    iteration_start_time = tic;
    % parfor
    for k = 1:Nt_coarse
        u_temp(:,k+1) = fine_solver(u_para(:,k),A,dx,dt_fine,J);
    end
    
    for k = 1:Nt_coarse
        u_para_new(:, k+1)=...
            coarse_solver(u_para_new(:,k),A,dx,dt_coarse)...
            - coarse_solver(u_para(:,k),A,dx,dt_coarse)...
            + u_temp(:,k+1);
    end
    u_para = u_para_new;
    
    u_err_para_seq = u_para(:,2:end) - u_seq(:,2:end);
    this_err = zeros(size(1:Nt_coarse));
    for k2 = 1 : Nt_coarse
        this_err(k2) = sqrt(1/Nx * u_err_para_seq(:,k2)' * u_err_para_seq(:,k2));
    end
    
    err(m+1) = max(this_err);
    
    % Assuming you already have x, u_para, T, Nt_coarse, m, and err(m) defined.

    % Flipping u_para by rows
    u_para_flipped = flipud(u_para);

    % Plot the figure
    % h_surf = surf(linspace(0, T, Nt_coarse+1), x, u_para_flipped, 'EdgeColor', 'none');
    % xlabel('Time (t)');
    % ylabel('Position (x)');
    % zlabel('u(x,t)');
    % title(['u(x,t) over time and position, m = ' num2str(m) ', err = ' num2str(err(m))]);
    % view(3); % 3D view

    % Set colormap and increase number of levels
    % colormap(jet(256)); % 'jet' colormap with 256 levels
    % colorbar;

    % Yticks_locs = get(gca, 'YTick');

    % Yticks_labels = L - Yticks_locs;

    % set(gca, 'YTickLabel', arrayfun(@num2str, Yticks_labels, 'UniformOutput', false));

    % pause(0.05);

    % clf;  % Use clf if you want to clear the figure for the next plot. Uncomment if needed.
    fprintf("Parallel iteration %d took %f seconds, %f minutes more \n", m, toc(iteration_start_time),toc(iteration_start_time)*(iter_max-m)/60);

end

% close(h);

figure(1);
m_list = 0:iter_max; 
semilogy(m_list, err, '-o');
hold on;
xlim([0,iter_max])

alpha1=zeros(iter_max-1, 1);
for i=1:1:iter_max-1
    alpha1(i)=exp((log(err(i+1))-log(err(i))));
end

%%
% Save the 'err' variable to a file in the 'plot' folder
% save('new J10/J_10_Radau_IIA_3_Solver.mat', 'err');
% save('new J10/J_10_Lobatto_IIIC_4_Solver.mat', 'err');
save('new J10/J_10_Lql_C4_Solver.mat', 'err');

% saveas(gcf, 'J20fig/J_20_Radau_IIA_3_Solver.png'); 
% saveas(gcf, 'J20fig/J_20_Lobatto_IIIC_4_Solver.png');
% saveas(gcf, 'J20fig/J_20_Lql_CR_Solver.png');

% Save the 'err' variable to a file in the 'plot' folder
% save('J100plot/J_100_Radau_IIA_3_Solver.mat', 'err');
% save('J100plot/J_100_Lobatto_IIIC_4_Solver.mat', 'err');
% save('J100plot/J_100_Lql_CR_Solver.mat', 'err');

% Save the 'err' variable to a file in the 'plot' folder
% save('v=0.02J2plot/J_2_Radau_IIA_3_Solver.mat', 'err');
% save('v=0.02J2plot/J_2_Lobatto_IIIC_4_Solver.mat', 'err');
% save('v=0.02J2plot/J_2_Lql_CR_Solver.mat', 'err');

% saveas(gcf, 'v=0.02J2fig/J_2_Radau_IIA_3_Solver.png'); 
% saveas(gcf, 'v=0.02J2fig/J_2_Lobatto_IIIC_4_Solver.png');
% saveas(gcf, 'v=0.02J2fig/J_2_Lql_CR_Solver.png');

% close(gcf);


disp(alpha1);


