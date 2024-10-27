function [H, b] = ComRPISet(H_phi, E, H_del, b_del, A_e)

k_max = 50;

H_r0 = H_phi*E;

k           = 1;
H_r_last    = H_r0;

while k < k_max
    k
    H_rk = [H_r_last; H_phi*E*A_e^k];

    H_r = H_rk;

    [nr, m] = size(H_r);

    opti = casadi.Opti();
    c    = opti.variable(nr, 1);
    d    = opti.variable(nr, 1);
    ksi  = opti.variable(m, 1);
    omgi = opti.variable(m, 1);

     opti.subject_to(H_r*ksi <= c + d);
     opti.subject_to(H_del*omgi <= b_del);
    for i = 1:nr
        opti.subject_to(c(i) <= H_r(i, :)*A_e*ksi)
        opti.subject_to(d(i) <= H_r(i, :)*omgi);
    end
    J = sum(c) + sum(d);
    opti.minimize(J);
    options = struct;
    options.ipopt.linear_solver = 'mumps';
    options.ipopt.print_level = 0;
    opti.solver('ipopt', options); 
    sol = opti.solve(); 

    % Extracting solutions
    c_sol = sol.value(c);
    d_sol = sol.value(d);
    if sol.value(J) == 0
        disp('Optimization successful');
        break
    else
        sum(c_sol) + sum(d_sol)
        k = k + 1;
        H_r_last = H_rk;
    end

end

H = H_rk;
b = c_sol + d_sol;

end




