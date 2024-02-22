function [maxh, tau_h] = compute_tau_h(parameters)
W          = parameters.W;
Psi        = parameters.Psi;
F_bar      = parameters.F_bar;
N          = parameters.N;
flag       = true;
n          = size(Psi,2);
m          = size(F_bar,1);
hh         = parameters.hh;
e_max      = parameters.e_max;
ogema_max  = W.omega_max;
gamma      = parameters.gamma;
alpha_max  = parameters.alpha_max;
tau_h      = N;

while flag
    A = F_bar*Psi^(tau_h + 1);
    for j = 1:m
        yalmip('clear')
        z = sdpvar(n,1);
        alpha = sdpvar(N+1,1); % a0, a1, ..., aN, totally N + 1
        cns=[];
        for k=1:N
            cns = [cns,F_bar*Psi^(k-1)*z <= ones(m,1)- alpha(k)*hh];
            cns = [cns,alpha(k)*e_max + ogema_max <= alpha(k+1)*ones(length(ogema_max),1)];
        end
        for k = N+1:tau_h+1
            cns = [cns,F_bar*Psi^(k-1)*z <= ones(m,1)-(gamma^(k-N-1)*(alpha(N+1) - alpha_max) + alpha_max)*hh];
        end
        obj = -A(j,:)*z - ((gamma^(tau_h+1-N)*(alpha(N+1) - alpha_max) + alpha_max)*hh(j));
        ops = sdpsettings('verbose',0);
        optimize(cns,obj,ops);
        maxh(j,1) = -value(obj);
    end
    if maxh <= ones(m,1)
        flag = false;
    else
       tau_h = tau_h + 1;
    end
end
end

