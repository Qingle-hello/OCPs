function unew = Lobatto_IIIC_3_Solver_BI(u,A,dx,tau,iter)
    for i=1:iter
        unew = Lobatto_IIIC_3_Solver_B(u,A,dx,tau);
        u = unew;
    end
end
