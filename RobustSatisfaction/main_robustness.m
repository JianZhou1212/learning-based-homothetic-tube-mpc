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

% Initial feasible region for LB Homothetic MPC
FeasibleSet = Homotheticfeasible_region(parameters, w_max_0);
Result_ht.FeasibleSet = FeasibleSet;

close all
figure(5)
plot(Result_ht.FeasibleSet, 'wire', 1, 'edgecolor', 'b', 'linewidth', 2);
title('Feasible Region')
legend('LB Homothetic', 'LB Rigid', 'Traditional Homothetic');
%%
V = Result_ht.FeasibleSet.V;
clc
N_MC = 30;
Feasible_Index = zeros(1, N_MC);

for h = 1:1:N_MC
    minValue = min(V(:, 1));
    maxValue = max(V(:, 2));

    minValueFormatted = round(minValue, 4);
    maxValueFormatted = round(maxValue, 4);
    
    state_feasible = [minValueFormatted; maxValueFormatted];

    Result_ht.x(:,1) = state_feasible;
    
    Iteration  = 30;
    new_sample = zeros(parameters.nx, 1);

    for i = 1:Iteration
        [Result_ht.s(:, i), Result_ht.c(:, i), Result_ht.alpha(:, i)] = HTMPC.recursive_mpc(Result_ht.x(:,i),Result_ht.w_max(:,i));
        [Result_ht.theta_k(:,i+1),Result_ht.rho_k(i+1),Result_ht.v_k(:,i+1), Result_ht.b] = HTMPC.update_W_hat(Result_ht.theta_k(:,i), Result_ht.rho_k(i), Result_ht.v_k(:,i), new_sample);
        Result_ht.u(:, i) = parameters.K*Result_ht.x(:,i) + Result_ht.c(1:parameters.nu,i);
        Result_ht.x(:,i + 1) = parameters.Phi*Result_ht.x(:,i) + parameters.B*Result_ht.c(1:parameters.nu,i) + new_sample;
        Result_ht.W_hat{i + 1} = Polyhedron(Result_ht.A, Result_ht.b);
        Result_ht.w_max(:,i+1) = Polysupport(Result_ht, (S.A)');
        
        new_sample = parameters.WSample_online(:, i);
    
    end

    X_Error = Result_ht.x(1, :);
    U_Input = Result_ht.u;
    exceed_up_control  = U_Input(U_Input > 3);
    exceed_low_control = U_Input(U_Input < -3);
    exceed_up_x        = X_Error(X_Error > 6.5);
    exceed_low_x       = X_Error(X_Error < -6.5);
    index = length(exceed_up_control) + length(exceed_low_control) + length(exceed_up_x) + length(exceed_low_x);
    if index == 0
        Feasible_Index(h) = 1;
    end

end
% Feasible Rate
sum(Feasible_Index)/N_MC
