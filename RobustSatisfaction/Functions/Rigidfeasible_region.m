function FeasibleSet=Rigidfeasible_region(parameters, Sk, hk, tau_k)
F_bar = parameters.F_bar;
Psi   = parameters.Psi;
nx    = parameters.nx;
nu    = parameters.nu;
N     = parameters.N;
cns   = [];
yalmip('clear')
Input_sample = sdpvar(nx,1);
xini  = sdpvar(nx,1);
s     = sdpvar(nx,1); 
c     = sdpvar(nu*N,1);
cns   = [cns, Sk.A*(xini - s) <= Sk.b];
for i = 1:tau_k+1
    cns = [cns,F_bar*Psi^(i-1)*[s;c]<= 1 - hk];
end
obj  = norm(xini - Input_sample);

ops  = sdpsettings('verbose',0);
MRPI = optimizer(cns, obj, ops, Input_sample, xini);

Sam_num = 3000;
Sample  = 50*rand(parameters.nx, Sam_num)-25;

parfor k = 1:Sam_num
     [sample_proj(:,k), ~] = MRPI(Sample(:,k));
end
Fs = sample_proj';
FeasibleSet = unique(Fs,'rows');
FeasibleSet = Polyhedron(FeasibleSet);

end
