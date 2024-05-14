function [uend, time1, cut, energy] = Lql_C2_Solver_AC(u0, Mass, Stiff, f, f_prime, Tend, tN, F, al, L_fast, U_fast)


    xN = length(Mass);
    cut = ones(tN, 1);
    energy = ones(tN+1,1);
    
    tau = Tend/tN;
    get_energy = @(u) (u' * Stiff * u) + ones(1, xN) * Mass * F(u);
    
    energy(1) = get_energy(u0);
    u(:,1) = u0;
    
    %if nargin == 9
        Lap_m = Mass \ Stiff;
    
        L = tau*Lap_m;
        I = diag(ones(size(u0)));   
        I_Euler_op = (0.37931*L^2 + 0.79033*L + I);      
        [L_IE,U_IE] = lu(I_Euler_op);
    %elseif nargin == 11
    %    L_IE = L_fast; U_IE = U_fast;
    %end
    
    t1=tic;
    for i=1:tN 
        
        f1 = f(u(:,i));
        
        r = (0.37447*L + I) * f1;
        b = (I-0.20967*L + 0.00484*L*L)* u(:,i) + tau * r;
        u(:,i+1) = (U_IE\(L_IE\(b)));   
        
        cut(i) = max(max(abs(min(max(u(:,i+1),-al),al) - u(:,i+1)))); 
        u(:,i+1) = min(max(u(:,i+1),-al),al);  
        energy(i+1) = get_energy(u(:,i+1));    
    
    end    
    time1 = toc(t1);
    
    uend = u(:,end);
    
      
    end
    