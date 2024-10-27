clear
close all
clc
addpath("Functions");

load('parameters2V.mat');
parameters = parameters2V;

W     = parameters.W;
Wtrue = parameters.Wtrue;
Phi   = parameters.Phi;
F     = parameters.F;
G     = parameters.G;
K     = parameters.K;
N     = parameters.N;
S     = parameters.S;
h     = parameters.h;
e_max = parameters.e_max;
w_max = parameters.w_max;
parameters.tau_h = N;

% initialization for W_hat of LB homothetic tube
HTMPC = HomotheticMPC(parameters);
[Result_ht.theta_k(:,1), Result_ht.rho_k(1), Result_ht.v_k(:,1)] = HTMPC.initialization_mpc( );
Result_ht.A = W.A;
Result_ht.b = Result_ht.theta_k(:,1) + (1-Result_ht.rho_k(1))*Result_ht.A*Result_ht.v_k(:,1);
Result_ht.W_hat{1}   = Polyhedron(Result_ht.A, Result_ht.b);
Result_ht.w_max(:,1) = Polysupport(Result_ht,(S.A)');
w_max_0              = Result_ht.w_max(:,1);

% compute horizon for LB homothetic MPC
[max_h, tau_h]       = compute_tau_h(parameters, w_max_0);
parameters.tau_h     = tau_h;

% compute horizon for traditional homothetic MPC
[max_th, tau_th]       = compute_tau_h(parameters, w_max);
parameters.tau_th     = tau_th;

% compute horizon for LB rigid tube
if max(h) < 1
tau_r = compute_tau_r(parameters, h);
Ps    = compute_Ps(parameters, h, tau_r);
end
parameters.tau_r = tau_r;
parameters.h     = h;
parameters.S     = S;
parameters.Ps    = Ps;
parameters.num_half_space_S = length(S.b);

% Initiazation for W_hat of LB RigidTube MPC
parameters.tau_k = N;
RTMPC = RigidMPC(parameters);
[Result_rt.rho_k(1), Result_rt.v_k(:, 1)] = RTMPC.initialization_mpc( );
Result_rt.A = W.A;
Result_rt.b = Result_rt.rho_k(:,1) + (1-Result_rt.rho_k(1))*Result_rt.A*Result_rt.v_k(:,1);
Result_rt.W_hat{1} = Polyhedron(Result_rt.A, Result_rt.b);

% Initial feasible region for LB Homothetic MPC
FeasibleSet = Homotheticfeasible_region(parameters, w_max_0);
Result_ht.FeasibleSet = FeasibleSet;

% Initial feasible region for Traditional Homothetic MPC
FeasibleSet_Tradition = TradHomotheticfeasible_region(parameters, w_max);
Result_th.FeasibleSet = FeasibleSet_Tradition;

% Initial feasible region for LB Ridig MPC
Sk    = Result_rt.rho_k(1)*S + (1 - Result_rt.rho_k(1))*(inv(eye(parameters.nx) - Phi)*Result_rt.v_k(:, 1));
hk    = Polysupport(Sk, (F + G*K)');
tau_k = compute_tau_k(parameters, hk);
Result_rt.FeasibleSet = Rigidfeasible_region(parameters, Sk, hk, tau_k);
%
close all
figure(5)
plot(Result_ht.FeasibleSet, 'wire', 1, 'edgecolor', 'b', 'linewidth', 2);
hold on
plot(Result_rt.FeasibleSet, 'wire', 1, 'edgecolor', 'k', 'linewidth', 2);
hold on
plot(Result_th.FeasibleSet, 'wire', 1, 'edgecolor', 'm', 'linewidth', 2);
title('Feasible Region')
legend('LB Homothetic', 'LB Rigid', 'Traditional Homothetic');

figure(6)
plot(Polyhedron(Result_ht.A, Result_ht.b), 'wire', 1, 'edgecolor', 'b', 'linewidth', 2);
hold on
plot(Polyhedron(Result_rt.A, Result_rt.b), 'wire', 1, 'edgecolor', 'k', 'linewidth', 2);
hold on
plot(Wtrue, 'wire', 1, 'edgecolor', 'm', 'linewidth', 2);
hold on
plot(parameters.W, 'wire', 1, 'edgecolor', 'r', 'linewidth', 2)
title('Control Set')
legend('Homothetic Tube MPC', 'Rigid Tube MPC', 'True', 'Robust');

figure(7)
plot(S, 'wire', 1, 'edgecolor', 'b', 'linewidth', 2);
hold on
plot(Sk, 'wire', 1, 'edgecolor', 'k', 'linewidth', 2);
title('RPMI Set')
legend('Original', 'Rigid Tube');
%%
clc
THMPC = TradHomotheticMPC(parameters);
state_feasible = [-6; 5];
state_infeasible = [-4; 6];
Result_ht.x(:,1) = state_feasible;
Result_th.x(:,1) = state_feasible;
Result_rt.x(:,1) = state_feasible;

Iteration  = 30;
new_sample = zeros(parameters.nx, 1);

SolTimeHT = ones(1, Iteration);
SolTimeRT = ones(1, Iteration);
SolTimeTH = ones(1, Iteration);

for i=1:Iteration

t0 = cputime;
[Result_ht.s(:, i), Result_ht.c(:, i), Result_ht.alpha(:, i)] = HTMPC.recursive_mpc(Result_ht.x(:,i),Result_ht.w_max(:,i));
[Result_ht.theta_k(:,i+1),Result_ht.rho_k(i+1),Result_ht.v_k(:,i+1), Result_ht.b] = HTMPC.update_W_hat(Result_ht.theta_k(:,i), Result_ht.rho_k(i), Result_ht.v_k(:,i), new_sample);
Result_ht.u(:, i) = parameters.K*Result_ht.x(:,i) + Result_ht.c(1:parameters.nu,i);
Result_ht.x(:,i + 1) = parameters.Phi*Result_ht.x(:,i) + parameters.B*Result_ht.c(1:parameters.nu,i) + new_sample;
Result_ht.W_hat{i + 1} = Polyhedron(Result_ht.A, Result_ht.b);
Result_ht.w_max(:,i+1) = Polysupport(Result_ht, (S.A)');
t1 = cputime-t0;
SolTimeHT(i) = t1;

t0 = cputime;
[Result_th.s(:, i), Result_th.c(:, i), Result_th.alpha(:, i)] = THMPC.recursive_mpc(Result_th.x(:,i));
Result_th.u(:, i)    = parameters.K*Result_th.x(:,i) + Result_th.c(1:parameters.nu,i);
Result_th.x(:,i + 1) = parameters.Phi*Result_th.x(:,i) + parameters.B*Result_th.c(1:parameters.nu,i) + new_sample;
t1 = cputime - t0;
SolTimeTH(i) = t1;


t0 = cputime;
Sk    = Result_rt.rho_k(i)*S + (1 - Result_rt.rho_k(i))*(inv(eye(parameters.nx) - Phi)*Result_rt.v_k(:, i));
hk    = Polysupport(Sk, (F + G*K)');
tau_k = compute_tau_k(parameters, hk);
parameters.tau_k = tau_k;
RTMPC = RigidMPC(parameters);
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
fprintf('Average solution time for THMPC is %.2f sec.\n', mean(SolTimeTH));
fprintf('Average solution time for RTMPC is %.2f sec.\n', mean(SolTimeRT));
%%
figure(1)
plot(Result_ht.x(1,:), Result_ht.x(2,:),'b-*', 'linewidth', 2.5)
hold on
plot(Result_rt.x(1,:), Result_rt.x(2,:),'k-p', 'linewidth', 2.5)
hold on
plot(Result_th.x(1,:), Result_th.x(2,:),'m-d', 'linewidth', 2.5)
title('State of Relative Model')
legend('Homothetic MPC', 'Rigid MPC', 'Traditional Homothetic MPC');

figure(2)
plot(Time, Result_ht.u(1, :),'b-*', 'linewidth', 2.5)
hold on
plot(Time, Result_rt.u(1, :),'k-p', 'linewidth', 2.5)
hold on
plot(Time, Result_th.u(1, :),'m-d', 'linewidth', 2.5)
title('Control Input of Relative Model')
legend('Homothetic MPC', 'Rigid MPC', 'Traditional Homothetic MPC');
%%
Result_Infeasible = struct;
Result_Infeasible.Result_ht = Result_ht;
Result_Infeasible.Result_th = Result_th;
Result_Infeasible.Result_rt = Result_rt;
Result_Infeasible.parameters = parameters;
save('Result_Infeasible.mat', 'Result_Infeasible');
%%
Result_Feasible = struct;
Result_Feasible.Result_ht = Result_ht;
Result_Feasible.Result_th = Result_th;
Result_Feasible.Result_rt = Result_rt;
Result_Feasible.parameters = parameters;
save('Result_Feasible.mat', 'Result_Feasible');



