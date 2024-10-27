function FeasibleSet=TradHomotheticfeasible_region(parameters, w_max)
yalmip('clear')
cns          = [ ];
Input_sample = sdpvar(parameters.nx, 1);
s            = sdpvar(parameters.nx, 1);
xini         = sdpvar(parameters.nx, 1);
c            = sdpvar(parameters.nu*parameters.N, 1);
alpha        = sdpvar(parameters.tau_th + 1, 1);
cns = [cns, alpha >= 0];
cns = [cns,parameters.S.A*(xini - s) <= alpha(1)*parameters.S.b];
for i   = 1:parameters.N
    cns = [cns, alpha(i)*parameters.e_max + w_max <= alpha(i+1)*ones(length(w_max),1)];
end
for i   = 1:parameters.tau_th + 1
    cns =[cns, parameters.F_bar*parameters.Psi^(i-1)*[s;c] <= 1 - alpha(i)*parameters.h];
end
for i   = parameters.N + 1:parameters.tau_th + 1
    cns =[cns, alpha(i) == 1];
end

obj  = norm(xini - Input_sample);
ops  = sdpsettings('verbose',0);
MRPI = optimizer(cns,obj,ops,Input_sample,xini);

Sam_num = 3000;
Sample  = 50*rand(parameters.nx, Sam_num) - 25;

parfor k = 1:Sam_num
     [sample_proj(:,k), ~] = MRPI(Sample(:,k));
end
Fs = sample_proj';
FeasibleSet = unique(Fs,'rows');
FeasibleSet = Polyhedron(FeasibleSet);

end
