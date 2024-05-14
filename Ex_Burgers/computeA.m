function A = computeA(Nx, nu, dx)
    main_diag = 2*nu/dx^2 * ones(Nx+1, 1);
    off_diag = -nu/dx^2 * ones(Nx, 1);

    A = diag(main_diag) + diag(off_diag, 1) + diag(off_diag, -1);

    A(1,1) = 1;  % Boundary conditions
    A(Nx+1, Nx+1) = 0;
end
