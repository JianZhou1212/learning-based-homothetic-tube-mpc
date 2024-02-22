clear
close all
clc
addpath("Functions");
systemmodel2V

load('parameters2V.mat');
parameters = parameters2V;

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
if  gamma<1
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

% compute rigid tube
epsilon = 0.01;
Sr      = compute_mrpi_set(Phi, W, epsilon);
hr      = Polysupport(Sr,(F + G*K)');
if max(hr) < 1
    tau_r = compute_tau_r(parameters, hr);
    Ps    = compute_Ps(parameters, hr, tau_r);
end
parameters.tau_r = tau_r;
parameters.hr    = hr;
parameters.Sr    = Sr;
parameters.Ps    = Ps;
parameters.num_half_space_S = length(Sr.b);

% Initialization for W_hat
HTMPC = HomotheticUQMPC(parameters);
[Result_ht.theta_k(:,1), Result_ht.rho_k(1), Result_ht.v_k(:,1)] = HTMPC.initialization_mpc( );
Result_ht.A = W.A;
Result_ht.b = Result_ht.theta_k(:,1) + (1-Result_ht.rho_k(1))*Result_ht.A*Result_ht.v_k(:,1);
Result_ht.W_hat{1} = Polyhedron(Result_ht.A, Result_ht.b);
Result_ht.omega_max(:,1) = Polysupport(Result_ht,(Sh.A)');
Result_ht.alpha_max(1)   = max(Result_ht.omega_max(:,1))/(1 - parameters.gamma);

% Initiazation for W_hat of RigidTube MPC
parameters.tau_k = 0;
RTMPC = RigidUQMPC(parameters);
[Result_rt.rho_k(1), Result_rt.v_k(:, 1)] = RTMPC.initialization_mpc( );
Result_rt.A = W.A;
Result_rt.b = Result_rt.rho_k(:,1) + (1-Result_rt.rho_k(1))*Result_ht.A*Result_rt.v_k(:,1);
Result_rt.W_hat{1} = Polyhedron(Result_rt.A, Result_rt.b);

% Initial feasible regions
FeasibleSet      = Homotheticfeasible_region(parameters, Result_ht.omega_max(:,1), Result_ht.alpha_max(1));
Result_ht.FeasibleSetHTMPC = Polyhedron(FeasibleSet);

Sk    = Result_rt.rho_k(1)*Sr + (1 - Result_rt.rho_k(1))*(inv(eye(parameters.nx) - Phi)*Result_rt.v_k(:, 1));
hk    = Polysupport(Sk,(F + G*K)');
tau_k = Com_tau_k(parameters, hk);
Result_rt.FeasibleSetRTMPC    = Rigidfeasible_region(parameters, Sk, hk, tau_k);

figure(5)
plot(Result_ht.FeasibleSetHTMPC, 'wire', 1, 'edgecolor', 'b', 'linewidth', 2);
hold on
plot(Result_rt.FeasibleSetRTMPC, 'wire', 1, 'edgecolor', 'k', 'linewidth', 2);
title('Feasible Region')
legend('Homothetic Tube MPC', 'Rigid Tube MPC');
%%
close all
clc
Result_ht.x(:,1) = [-12; 6];
Result_rt.x(:,1) = [-12; 6];

Iteration  = 30;
new_sample = zeros(parameters.nx, 1);
SolTimeHT = ones(1, Iteration);
SolTimeRT = ones(1, Iteration);
for i=1:Iteration

    i

    t0 = cputime;
    [Result_ht.s(:, i), Result_ht.c(:, i), Result_ht.alpha(:, i)]= HTMPC.recursive_mpc(Result_ht.x(:,i),Result_ht.omega_max(:,i),Result_ht.alpha_max(i));
    [Result_ht.theta_k(:,i+1),Result_ht.rho_k(i+1),Result_ht.v_k(:,i+1), Result_ht.b, Result_ht.omega_max(:,i+1), Result_ht.alpha_max(i+1)] = HTMPC.update_W_hat(Result_ht.theta_k(:,i), Result_ht.rho_k(i), Result_ht.v_k(:,i), new_sample);
    Result_ht.u(:, i) = parameters.K*Result_ht.x(:,i) + Result_ht.c(1:parameters.nu,i);
    Result_ht.x(:,i + 1) = parameters.Phi*Result_ht.x(:,i) + parameters.B*Result_ht.c(1:parameters.nu,i) + new_sample;
    Result_ht.W_hat{i + 1}  = Polyhedron(Result_ht.A, Result_ht.b);
    t1 = cputime-t0;
    SolTimeHT(i) = t1;

    t0 = cputime;
    Sk    = Result_rt.rho_k(i)*Sr + (1 - Result_rt.rho_k(i))*(inv(eye(parameters.nx) - Phi)*Result_rt.v_k(:, i));
    hk    = Polysupport(Sk,(F + G*K)');
    tau_k = Com_tau_k(parameters, hk);
    parameters.tau_k = tau_k;
    Result_rt.Sk_time{i} = Sk;
    [Result_rt.s(:, i), Result_rt.c(:, i)] = RTMPC.recursive_mpc(Result_rt.x(:,i), Sk, hk);
    [Result_rt.rho_k(i+1), Result_rt.v_k(:, i+1)]= RTMPC.update_W_hat(Result_rt.rho_k(i), Result_rt.v_k(:, i), new_sample);
    Result_rt.b = Result_rt.rho_k(:,i+1) + (1-Result_rt.rho_k(i+1))*Result_ht.A*Result_rt.v_k(:,i+1);
    Result_rt.W_hat{i + 1} = Polyhedron(Result_rt.A, Result_rt.b);
    Result_rt.u(:, i) = parameters.K*Result_rt.x(:,i) + Result_rt.c(1:parameters.nu,i);
    Result_rt.x(:,i+1)  = parameters.Phi*Result_rt.x(:,i) + parameters.B*Result_rt.c(1:parameters.nu,i) + new_sample;
    t1 = cputime-t0;
    SolTimeRT(i) = t1;

    new_sample = parameters.WSample_online(:, i);

end
Time = 0:parameters.T:parameters.T*(Iteration-1);
fprintf('Average solution time for HTMPC is %.2f sec.\n', mean(SolTimeHT));
fprintf('Average solution time for RTMPC is %.2f sec.\n', mean(SolTimeRT));

figure(1)
plot(Result_ht.x(1,:), Result_ht.x(2,:),'b-*', 'linewidth', 2.5)
hold on
plot(Result_rt.x(1,:), Result_rt.x(2,:),'k-p', 'linewidth', 2.5)
title('State of Relative Model')
legend('Homothetic Tube MPC', 'Rigid Tube MPC');

figure(2)
plot(Time, Result_ht.u(1, :),'b-*', 'linewidth', 2.5)
hold on
plot(Time, Result_rt.u(1, :),'k-p', 'linewidth', 2.5)
title('Control Input of Relative Model')
legend('Homothetic Tube MPC', 'Rigid Tube MPC');

figure(3)
plot(Result_ht.alpha(1, 1)*Polyhedron(Sh.A, Sh.b), 'wire', 1, 'edgecolor', 'b', 'linewidth', 2.5);
hold on
plot(Result_rt.Sk_time{1}, 'wire', 1, 'edgecolor', 'k', 'linewidth', 2.5);
hold on
plot(Sr, 'wire', 1, 'edgecolor', 'r', 'linewidth', 2.5);
title('Tube')
legend('Homothetic Tube MPC', 'Rigid Tube MPC', 'MRPI Set');

figure(4)
plot(Result_ht.W_hat{1}, 'wire', 1, 'edgecolor', 'b', 'linewidth', 2.5);
hold on
plot(Result_rt.W_hat{1}, 'wire', 1, 'edgecolor', 'k', 'linewidth', 2.5);
hold on
plot(Polyhedron(parameters.W.A, parameters.W.b), 'wire', 1, 'edgecolor', 'r', 'linewidth', 2.5)
title('Disturbance Set');
legend('Homothetic Tube MPC', 'Rigid Tube MPC', 'Max. Disturbanc Set');

% Result_INFeasible = struct;
% Result_INFeasible.Result_ht = Result_ht;
% Result_INFeasible.Result_rt = Result_rt;
% Result_INFeasible.parameters = parameters;
% save('Result_INFeasible.mat', 'Result_INFeasible');
% 
Result_Feasible = struct;
Result_Feasible.Result_ht = Result_ht;
Result_Feasible.Result_rt = Result_rt;
Result_Feasible.parameters = parameters;
save('Result_Feasible.mat', 'Result_Feasible');



