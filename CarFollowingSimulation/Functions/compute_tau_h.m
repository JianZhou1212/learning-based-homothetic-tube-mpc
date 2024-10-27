function [maxh, n] = compute_tau_h(parameters, w_max)
% Compute tau by Lemma 1
Psi        = parameters.Psi;
F_bar      = parameters.F_bar;
N          = parameters.N;
flag       = true;
h          = parameters.h;
e_max      = parameters.e_max;
n          = N;

while flag
    A = F_bar*Psi^(n + 1);
    for j = 1:size(F_bar,1)
        opti  = casadi.Opti();
        z     = opti.variable(size(Psi,2), 1);
        alpha = opti.variable(n + 1, 1); % a0, a1, ..., an, totally n + 1
        for k = 1:N % 0,...,N-1
            opti.subject_to(alpha(k)*e_max + w_max <= alpha(k + 1)*ones(length(w_max), 1));
        end
        for k = 1:n + 1 % 0, ..., n
            opti.subject_to(F_bar*Psi^(k-1)*z <= ones(size(F_bar,1),1)- alpha(k)*h);
        end
        for k = N + 1:n + 1 % N, ..., n
            opti.subject_to(alpha(k) == 1);
        end
        opti.subject_to(alpha >= 0);
        J = -A(j,:)*z;
        opti.minimize(J);
        options = struct;
        options.ipopt.linear_solver = 'ma57';
        options.ipopt.print_level = 0;
        opti.solver('ipopt', options); 
        sol = opti.solve(); 

        maxh(j,1) = -sol.value(J);
    end
    if maxh <= 1 - h
        flag = false;
    else
       n = n + 1;
    end
end
end

