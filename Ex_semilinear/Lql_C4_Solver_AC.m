function [uend, time1, cut, energy] = Lql_C4_Solver_AC(u0, Mass, Stiff, f, f_prime, Tend, tN, F, al, L_fast, U_fast)

    [NoF,~] = size(u0);
    
    
    xN = length(Mass);
    cut = ones(tN, 1);
    energy = ones(tN+1,1);
    
    
    tau = Tend/tN;
    get_energy = @(u) (u' * Stiff * u) + ones(1, xN) * Mass * F(u);
    
    energy(1) = get_energy(u0);
    u(:,1) = u0;
    
    % if nargin == 9
        Lap_m = Mass \ Stiff;
    
        L = tau*Lap_m;
        I = diag(ones(size(u0)));   
        I_Euler_op = (0.37898*L^2 + 0.78909*L + I);  
        A = inv(I_Euler_op);
    % elseif nargin == 11
    %     L_IE = L_fast; U_IE = U_fast;
    % end
    
    t1 = tic;
    for i=1:tN 
        f1 = f(u(:,i));
        
        r = (0.37422*L + I) * f1;
        b = (0.98744*I - 0.22082*L) * u(:,i) + tau * r;
    
        u(:,i+1) = A * b +    0.01256 * u(:,i);
        
        cut(i) = max(max(abs(min(max(u(:,i+1),-al),al) - u(:,i+1)))); 
        u(:,i+1) = min(max(u(:,i+1),-al),al);  
        energy(i+1) = get_energy(u(:,i+1));    
    
    end    
    
    time1 = toc(t1);
    
    uend = u(:,end);
    
        
    
    end
    