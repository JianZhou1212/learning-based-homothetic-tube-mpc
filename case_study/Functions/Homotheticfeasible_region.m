function FeasibleSet=Homotheticfeasible_region(parameters, omega_max, alpha_max)
yalmip('clear')
cns          = [];
Input_sample = sdpvar(parameters.nx,1);
s            = sdpvar(parameters.nx,1);
xini         = sdpvar(parameters.nx,1);
c            = sdpvar(parameters.nu*parameters.N,1);
alpha        = sdpvar(parameters.N + 1,1);
cns = [cns, alpha >= 0];
cns = [cns,parameters.Sh.A*(xini - s) <= alpha(1)*parameters.Sh.b];
for i = 1:parameters.N
    cns = [cns,parameters.F_bar*parameters.Psi^(i-1)*[s; c] <= 1 - alpha(i)*parameters.hh];
    cns = [cns,alpha(i)*parameters.e_max + omega_max <= alpha(i+1)*ones(length(omega_max),1)];
end

for i = parameters.N + 1:parameters.tau_h + 1
    cns =[cns,parameters.F_bar*parameters.Psi^(i-1)*[s;c] <= 1-(parameters.gamma^(i-parameters.N-1)*(alpha(parameters.N+1)- alpha_max) + alpha_max)*parameters.hh];
end

obj  = norm(xini - Input_sample);
ops  = sdpsettings('verbose',0);
MRPI = optimizer(cns,obj,ops,Input_sample,xini);

Sam_num = 5000;
Sample  = 50*rand(parameters.nx,Sam_num) - 25;

parfor kkk = 1:Sam_num
     [sample_proj(:,kkk),errorcode] = MRPI(Sample(:,kkk));
end
Fs = sample_proj';
FeasibleSet = unique(Fs,'rows');

end
