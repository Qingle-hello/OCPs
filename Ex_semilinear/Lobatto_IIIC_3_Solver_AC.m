function [uend, time1, cut, energy] = Lobatto_IIIC_3_Solver_AC(u0, Mass, Stiff, f, f_prime, Tend, tN, F, al)

    xN = length(Mass);
    cut = ones(tN, 1);
    energy = ones(tN+1,1);
    
    
    tau = Tend/tN;
    get_energy = @(u) (u' * Stiff * u) + ones(1, xN) * Mass * F(u);
    
    energy(1) = get_energy(u0);
    u(:,1) = u0;
    
    Lap_m = Mass \ Stiff;
    
    I = diag(ones(size(u0)));
    
    RHS = @(u) - Lap_m * u + f(u);
    RHS_prime = @(u) - Lap_m + diag(f_prime(u));

    a11 = 1/6;a12=-1/3;a13=1/6;
    a21=1/6;a22=5/12;a23=-1/12;
    a31=1/6;a32=2/3;a33=1/6;
    b1=1/6;b2=2/3;b3=1/6;

    t1=tic;
    for i= 1:tN

        ques = @(uu) [ ...
            u(:,i) + tau*a11*RHS(uu(1:end/3)) + tau*a12*RHS(uu((end/3+1):(2*end/3))) + tau*a13*RHS(uu((2*end/3+1):end)) - uu(1:end/3); ...
            u(:,i) + tau*a21*RHS(uu(1:end/3)) + tau*a22*RHS(uu((end/3+1):(2*end/3))) + tau*a23*RHS(uu((2*end/3+1):end)) - uu((end/3+1):(2*end/3)); ...
            u(:,i) + tau*a31*RHS(uu(1:end/3)) + tau*a32*RHS(uu((end/3+1):(2*end/3))) + tau*a33*RHS(uu((2*end/3+1):end)) - uu((2*end/3+1):end)] ;
        ques_prime = @(uu) [ ...
            tau*a11*RHS_prime(uu(1:end/3)) - I, tau*a12*RHS_prime(uu((end/3+1):(2*end/3))), tau*a13*RHS_prime(uu((2*end/3+1):end));...
            tau*a21*RHS_prime(uu(1:end/3)), tau*a22*RHS_prime(uu((end/3+1):(2*end/3))) - I, tau*a23*RHS_prime(uu((2*end/3+1):end));...
            tau*a31*RHS_prime(uu(1:end/3)), tau*a32*RHS_prime(uu((end/3+1):(2*end/3))), tau*a33*RHS_prime(uu((2*end/3+1):end)) - I];
        x0 = [u(:,i);u(:,i);u(:,i)];
        
        answer = Newton_iter(ques, ques_prime, x0, 1e-12);
        
        u(:,i+1) = u(:,i) + tau*b1*RHS(answer(1:end/3)) + tau*b2*RHS(answer((end/3+1):(2*end/3))) + tau*b3*RHS(answer((2*end/3+1):end));
    %     u(:,i+1) = fsolve(ques, u(:,i), opt);
    %     u(:,i+1) = r1;
         
        u(:,i+1) = min(max(u(:,i+1),-al),al);  

    end
    time1=toc(t1);
    
    
    uend = u(:,end);
    end