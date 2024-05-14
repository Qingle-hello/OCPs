function [Mass, Stiff, x, NoF] = fem_matrix_nm(L, M, r)
    % r refers to the order of FEM
    % L refers to the length of the interval
    % M refers to the count of the subintervals, or mesh sizes.
    %FEM_MATRIX Summary of this function goes here
    %   Detailed explanation goes here
    NoF = M * r + 1; x = zeros(NoF, 1); h = L / M;
    [mass_local, stiff_local, xg, wg] = local_matrix(L, M, r);
    Mass = sparse(NoF, NoF); Stiff = sparse(NoF, NoF);
    %x(1) = -L / 2;

    for k = 1:M
        Mass(((k - 1) * r + 1):(k * r + 1), ((k - 1) * r + 1):(k * r + 1)) = Mass(((k - 1) * r + 1):(k * r + 1), ((k - 1) * r + 1):(k * r + 1)) + mass_local;
        Stiff(((k - 1) * r + 1):(k * r + 1), ((k - 1) * r + 1):(k * r + 1)) = Stiff(((k - 1) * r + 1):(k * r + 1), ((k - 1) * r + 1):(k * r + 1)) + stiff_local;
        x(((k - 1) * r + 2):(k * r + 1)) = h * (k - 1) + (xg(2:end)' + 1) * h / 2;
    end

    % Mass(1,:)=[];Stiff(1,:)=[];
    % Mass(:,1)=[];Stiff(:,1)=[];
    % Mass(end,:)=[];Stiff(end,:)=[];
    % Mass(:,end)=[];Stiff(:,end)=[];
    % NoF=NoF-2;

end
