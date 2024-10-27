function tau_k = compute_tau_k(parameters, hk) 
    zeta = max((1 - hk)./(1 - parameters.h));
    tau_k = parameters.tau_r;
    vector = [ ];
    for i = 1:1:parameters.nc
        vector = [vector parameters.F_bar(i, :)*parameters.Psi^(tau_k + 1)*pinv(parameters.Ps)*(parameters.Psi^(tau_k + 1))'*(parameters.F_bar(i, :))' - (1 - hk(i))/zeta^2];
    end
    while max(vector) > 0
        tau_k = tau_k + 1;
        vector = [ ];
        for i = 1:1:parameters.nc
            vector = [vector parameters.F_bar(i, :)*parameters.Psi^(tau_k + 1)*pinv(parameters.Ps)*parameters.Psi^(tau_k + 1)*(parameters.F_bar(i, :))' - (1 - hk(i))/zeta^2];
        end
    end
end