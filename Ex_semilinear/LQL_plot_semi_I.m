close all; clear all;clc;


figure(1); %hold on;
%%
subplot(2,2,1);
Lb3c2N2_20 = load("Ex_2_a_data/figure data/J20/J_20_Lobatto_IIIC_2_Solver.mat").err;
Lql_solver_20 = load("Ex_2_a_data/figure data/J20/J_20_Lql_C2_Solver.mat").err;

Lb3c2N2_100 = load("Ex_2_a_data/figure data/J100/J_100_Lobatto_IIIC_2_Solver.mat").err;
Lql_solver_100 = load("Ex_2_a_data/figure data/J100/J_100_Lql_C2_Solver.mat").err;

Lb3c2N2_20 = Lb3c2N2_20(1:11);
Lql_solver_20 = Lql_solver_20(1:11);
Lb3c2N2_100 = Lb3c2N2_100(1:11);
Lql_solver_100 = Lql_solver_100(1:11);


m_list = 0:10;
semilogy(m_list, Lb3c2N2_20, '-^');hold on;
semilogy(m_list, Lql_solver_20, '-o');
semilogy(m_list, Lb3c2N2_100, '-^');
semilogy(m_list, Lql_solver_100, '-o');

xlim([0,10])
ref = zeros(size(m_list));
for i = m_list+1
    ref(i) = (0.3)^((i-1));
end
semilogy(m_list(1:5), ref(1:5), '--');

ref1 = zeros(size(m_list));
for i = m_list+1
    ref1(i) = (0.04)^((i-1));
end
semilogy(m_list(1:5), ref1(1:5), '--');
ylim([1e-15,1])

legend("J=20 BE","J=20 optimized","J=100 BE","J=100 optimized","ref 0.3","ref 0.04", 'Location', 'southwest');
title("Lobatto IIIC, 2 stage");
xlabel("iteration");
ylabel("error");


%%
subplot(2,2,2);
Lb3c2N2_20 = load("Ex_2_a_data/figure data/J20/J_20_Lobatto_IIIC_3_Solver.mat").err;
Lql_solver_20 = load("Ex_2_a_data/figure data/J20/J_20_Lql_C3_Solver.mat").err;

Lb3c2N2_100 = load("Ex_2_a_data/figure data/J100/J_100_Lobatto_IIIC_3_Solver.mat").err;
Lql_solver_100 = load("Ex_2_a_data/figure data/J100/J_100_Lql_C3_Solver.mat").err;


Lb3c2N2_20 = Lb3c2N2_20(1:11);
Lql_solver_20 = Lql_solver_20(1:11);
Lb3c2N2_100 = Lb3c2N2_100(1:11);
Lql_solver_100 = Lql_solver_100(1:11);


m_list = 0:10;
semilogy(m_list, Lb3c2N2_20, '-^');hold on;
semilogy(m_list, Lql_solver_20, '-o');
semilogy(m_list, Lb3c2N2_100, '-^');
semilogy(m_list, Lql_solver_100, '-o');

xlim([0,10])
ref = zeros(size(m_list));
for i = m_list+1
    ref(i) = (0.3)^((i-1));
end
semilogy(m_list(1:5), ref(1:5), '--');

ref1 = zeros(size(m_list));
for i = m_list+1
    ref1(i) = (0.04)^((i-1));
end
semilogy(m_list(1:5), ref1(1:5), '--');
ylim([1e-15,1])

legend("J=20 BE","J=20 optimized","J=100 BE","J=100 optimized","ref 0.3","ref 0.04", 'Location', 'southwest');
title("Lobatto IIIC, 3 stage");
xlabel("iteration");
ylabel("error");

%%
subplot(2,2,3);
Lb3c2N2_20 = load("Ex_2_a_data/figure data/J20/J_20_Lobatto_IIIC_4_Solver.mat").err;
Lql_solver_20 = load("Ex_2_a_data/figure data/J20/J_20_Lql_C4_Solver.mat").err;

Lb3c2N2_100 = load("Ex_2_a_data/figure data/J100/J_100_Lobatto_IIIC_4_Solver.mat").err;
Lql_solver_100 = load("Ex_2_a_data/figure data/J100/J_100_Lql_C4_Solver.mat").err;


Lb3c2N2_20 = Lb3c2N2_20(1:11);
Lql_solver_20 = Lql_solver_20(1:11);
Lb3c2N2_100 = Lb3c2N2_100(1:11);
Lql_solver_100 = Lql_solver_100(1:11);


m_list = 0:10;
semilogy(m_list, Lb3c2N2_20, '-^');hold on;
semilogy(m_list, Lql_solver_20, '-o');
semilogy(m_list, Lb3c2N2_100, '-^');
semilogy(m_list, Lql_solver_100, '-o');

xlim([0,10])
ref = zeros(size(m_list));
for i = m_list+1
    ref(i) = (0.3)^((i-1));
end
semilogy(m_list(1:5), ref(1:5), '--');

ref1 = zeros(size(m_list));
for i = m_list+1
    ref1(i) = (0.04)^((i-1));
end
semilogy(m_list(1:5), ref1(1:5), '--');
ylim([1e-15,1])

legend("J=20 BE","J=20 optimized","J=100 BE","J=100 optimized","ref 0.3","ref 0.04", 'Location', 'southwest');
title("Lobatto IIIC, 4 stage");
xlabel("iteration");
ylabel("error");

%%
subplot(2,2,4);
Lb3c2N2_20 = load("Ex_2_a_data/figure data/J20/J_20_Radau_IIA_3_Solver.mat").err;
Lql_solver_20 = load("Ex_2_a_data/figure data/J20/J_20_Lql_CR_Solver.mat").err;

Lb3c2N2_100 = load("Ex_2_a_data/figure data/J100/J_100_Radau_IIA_3_Solver.mat").err;
Lql_solver_100 = load("Ex_2_a_data/figure data/J100/J_100_Lql_CR_Solver.mat").err;


Lb3c2N2_20 = Lb3c2N2_20(1:11);
Lql_solver_20 = Lql_solver_20(1:11);
Lb3c2N2_100 = Lb3c2N2_100(1:11);
Lql_solver_100 = Lql_solver_100(1:11);


m_list = 0:10;
semilogy(m_list, Lb3c2N2_20, '-^');hold on;
semilogy(m_list, Lql_solver_20, '-o');
semilogy(m_list, Lb3c2N2_100, '-^');
semilogy(m_list, Lql_solver_100, '-o');

xlim([0,10])
ref = zeros(size(m_list));
for i = m_list+1
    ref(i) = (0.3)^((i-1));
end
semilogy(m_list(1:5), ref(1:5), '--');

ref1 = zeros(size(m_list));
for i = m_list+1
    ref1(i) = (0.04)^((i-1));
end
semilogy(m_list(1:5), ref1(1:5), '--');
ylim([1e-15,1])

legend("J=20 BE","J=20 optimized","J=100 BE","J=100 optimized","ref 0.3","ref 0.04", 'Location', 'southwest');
title("Radau IIA, 3 stage");
xlabel("iteration");
ylabel("error");
%%
set(1,'position',[100 100 900 600])
saveas(1, "fig103", "pdf")
