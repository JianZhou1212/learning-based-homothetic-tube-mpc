function n = compute_tau_r(system, h)
Psi   = system.Psi;
F_bar = system.F_bar;
N     = system.N;
flag  = true;
m     = size(F_bar,1);
n     = N;
while flag
    A = F_bar*Psi^(n+1);
    for j = 1:m
        opti = casadi.Opti();
        z    = opti.variable(size(Psi,2),1);
        for k   = 1:n+1
            opti.subject_to(F_bar*Psi^(k-1)*z <= ones(m,1) - h);
        end
        J = -A(j,:)*z;
        opti.minimize(J);
        options = struct;
        options.ipopt.linear_solver = 'ma57';
        options.ipopt.print_level = 0;
        opti.solver('ipopt', options); 
        sol       = opti.solve(); 
        maxh(j,1) = -sol.value(J);
    end
    if maxh <= 1 - h
        flag = false;
    else
        n = n+1;
        n
    end
end
end

