clear
close all
clc
addpath("Functions");
systemmodel5V

load('parameters5V.mat');
parameters = parameters5V;

W    = parameters.W;
Phi  = parameters.Phi;
Fx   = parameters.Fx;
fx   = parameters.fx;
F    = parameters.F;
G    = parameters.G;
K    = parameters.K;

% compute homothetic tube
[Sh.A, Sh.b] = findMRPIset(Phi, Fx, fx, W.A, W.b); 
e_max        = Polysupport(Sh,(Sh.A*Phi)'); 
gamma        = max(e_max);
hh           = Polysupport(Sh,(F + G*K)');     
W.omega_max  = Polysupport(W,(Sh.A)');  
W.mu         = max(W.omega_max);          
parameters.W = W;
if  gamma < 1
    parameters.Sh           = Sh;
    parameters.e_max        = e_max;
    parameters.gamma        = gamma;
    parameters.alpha_max    = W.mu/(1-gamma);
    parameters.hh           = hh; 
   [maxh, parameters.tau_h] = compute_tau_h(parameters); % maxh should be leq. 1
end

gamma = parameters.gamma;
N = parameters.N;
qa = 0.1;
pa = 0.1 + qa/(1 - gamma^2);
Pa = diag([qa*ones(1, N) pa]);
parameters.Pa = Pa; % Cost function (31)

% Initialization for W_hat
HTMPC = HomotheticUQMPC(parameters);
[Result_ht.theta_k(:,1), Result_ht.rho_k(1), Result_ht.v_k(:,1)] = HTMPC.initialization_mpc( );
Result_ht.A = W.A;
Result_ht.b = Result_ht.theta_k(:,1) + (1-Result_ht.rho_k(1))*Result_ht.A*Result_ht.v_k(:,1);
Result_ht.W_hat{1} = Polyhedron(Result_ht.A, Result_ht.b);
Result_ht.omega_max(:,1) = Polysupport(Result_ht,(Sh.A)');
Result_ht.alpha_max(1)   = max(Result_ht.omega_max(:,1))/(1 - parameters.gamma);
%%
% Initial feasible regions
FeasibleSet      = Homotheticfeasible_region(parameters, Result_ht.omega_max(:,1), Result_ht.alpha_max(1));
Result_ht.FeasibleSetHTMPC = Polyhedron(FeasibleSet);
%%
close all
clc
Result_ht.x(:,1) = [-12; 6; -12; 4; -12; 2; -12; 0];

Iteration  = 30;
new_sample = zeros(parameters.nx, 1);
SolTimeHT = ones(1, Iteration);
SolTimeRT = ones(1, Iteration);
for i=1:Iteration

    t0 = cputime;
    [Result_ht.s(:, i), Result_ht.c(:, i), Result_ht.alpha(:, i)]= HTMPC.recursive_mpc(Result_ht.x(:,i),Result_ht.omega_max(:,i),Result_ht.alpha_max(i));
    [Result_ht.theta_k(:,i+1),Result_ht.rho_k(i+1),Result_ht.v_k(:,i+1), Result_ht.b, Result_ht.omega_max(:,i+1), Result_ht.alpha_max(i+1)] = HTMPC.update_W_hat(Result_ht.theta_k(:,i), Result_ht.rho_k(i), Result_ht.v_k(:,i), new_sample);
    Result_ht.u(:, i) = parameters.K*Result_ht.x(:,i) + Result_ht.c(1:parameters.nu,i);
    Result_ht.x(:,i + 1) = parameters.Phi*Result_ht.x(:,i) + parameters.B*Result_ht.c(1:parameters.nu,i) + new_sample;
    Result_ht.W_hat{i + 1}  = Polyhedron(Result_ht.A, Result_ht.b);
    t1 = cputime-t0;
    SolTimeHT(i) = t1;

    new_sample = parameters.WSample_online(:, i);

end
Time = 0:parameters.T:parameters.T*(Iteration-1);
fprintf('Average solution time for HTMPC is %.2f sec.\n', mean(SolTimeHT));

close all
figure(1)
plot(Result_ht.x(1,:), Result_ht.x(2,:),'b-*', 'linewidth', 2.5)
hold on
plot(Result_ht.x(3,:), Result_ht.x(4,:),'k-p', 'linewidth', 2.5)
hold on
plot(Result_ht.x(5,:), Result_ht.x(6,:),'g-p', 'linewidth', 2.5)
hold on
plot(Result_ht.x(7,:), Result_ht.x(8,:),'r-*', 'linewidth', 2.5)
title('State of Relative Model')
legend('1--2 Vehicle', '2--3 Vehicle', '3--4 Vehicle', '4--5 Vehicle');

figure(2)
plot(Time, Result_ht.u(1, :),'b-*', 'linewidth', 2.5)
hold on
plot(Time, Result_ht.u(2, :),'k-p', 'linewidth', 2.5)
hold on
plot(Time, Result_ht.u(3, :),'g-p', 'linewidth', 2.5)
hold on
plot(Time, Result_ht.u(4, :),'r-*', 'linewidth', 2.5)
title('Control Input of Relative Model')
legend('Vehicle 2', 'Vehicle 3', 'Vehicle 4', 'Vehicle 5');

Result_Feasible5V = struct;
Result_Feasible5V.Result_ht = Result_ht;
Result_Feasible5V.parameters = parameters;
save('Result_Feasible5V.mat', 'Result_Feasible5V');



