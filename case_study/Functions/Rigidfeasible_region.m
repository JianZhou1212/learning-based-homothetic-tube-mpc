function FeasibleSet = Rigidfeasible_region(parameters, Sk, hk, tau_k)
    Sam_num = 5000;
    Sample = 50*rand(2, Sam_num) - 25;
    yalmip('clear')
    z = sdpvar(size(parameters.Psi, 2), 1);
    proj_sample =  sdpvar(parameters.nx, 1);
    Input_sample = sdpvar(parameters.nx, 1);
    cns = [];
    for i = 0:1:tau_k
        cns = [cns,parameters.F_bar*(parameters.Psi^i)*z <= ones(parameters.nc, 1) - hk];
    end
    cns = [cns,z(1:parameters.nx) == proj_sample];
    J = (proj_sample - Input_sample)'*(proj_sample - Input_sample);
    ops = sdpsettings('relax', 0);
    Fea_Set = optimizer(cns, J, ops, Input_sample, proj_sample);
    
    % evaluate the function and form a convexhull
    parfor k = 1:Sam_num
        [sample_proj(:,k),~] = Fea_Set(Sample(:,k));
    end
    [convhull_index, ~] = convhull(sample_proj');
    F_Ns_hat_opt = sample_proj(:,convhull_index);
    F_Ns_hat_opt = Polyhedron(F_Ns_hat_opt');
    FeasibleSet = F_Ns_hat_opt + Sk;
end