function Ps=compute_Ps(system, hr, taur)
Psi   = system.Psi;
F_bar = system.F_bar;

m     = size(F_bar,1);
nx    = system.nx;
nu    = system.nu;
N     = system.N;

Sam_num = 1000;
Sample  = 50*rand(nx+nu*N, Sam_num)-25;

yalmip('clear')
proj_sample  = sdpvar(nx + nu*N, 1);
Input_sample = sdpvar(nx + nu*N, 1);
cns = [ ];
for k   = 1:taur
    cns = [cns, F_bar*Psi^k*proj_sample <= ones(m, 1) - hr];
end
obj = norm(proj_sample - Input_sample);
ops = sdpsettings('verbose',0);
OPT = optimizer(cns, obj, ops, Input_sample, proj_sample);

parfor k = 1:Sam_num
     [sample_proj(:, k), errorcode] = OPT(Sample(:, k));
end
Fs = sample_proj';
Fs = unique(Fs,'rows');

yalmip('clear')
Ps      = sdpvar(nx + N*nu, nx + N*nu); 
cns     = [ ];
for i   = 1:1:size(Fs, 1)
    cns = [cns, norm(Ps*Fs(i, :)', 2) <= 1]; 
end

obj = -logdet(Ps); 
ops = sdpsettings('verbose',0);
sol = optimize(cns, obj, ops);
Ps  = value(Ps);
end

