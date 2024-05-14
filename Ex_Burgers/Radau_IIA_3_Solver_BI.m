function unew = Radau_IIA_3_Solver_BI(u,A,dx,tau,iter)
    for i=1:iter
        unew = Radau_IIA_3_Solver_B(u,A,dx,tau);
        u = unew;
    end
end
