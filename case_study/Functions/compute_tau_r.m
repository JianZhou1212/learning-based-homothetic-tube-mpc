function tau_r=compute_tau_r(system, hr)
Psi   = system.Psi;
F_bar = system.F_bar;
N     = system.N;
flag  = true;
n     = size(Psi,2);
m     = size(F_bar,1);
tau_r = N;
while flag
    A = F_bar*Psi^(tau_r+1);
    for j = 1:m
        yalmip('clear')
        z   = sdpvar(n,1);
        cns = [];
        for k   = 1:tau_r+1
            cns = [cns,F_bar*Psi^(k-1)*z <= ones(m,1) - hr];
        end
        obj = -A(j,:)*z;
        ops = sdpsettings('verbose',0);
        optimize(cns,obj,ops);
        maxh(j,1) = A(j,:)*value(z);
    end
    if maxh <= ones(m,1)-hr
        flag = false;
    else
        tau_r = tau_r+1;
    end
end
end

