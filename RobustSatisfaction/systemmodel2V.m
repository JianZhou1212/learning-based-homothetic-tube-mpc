clc
clear 
close all
addpath("Functions");

% system parameters
T  = 1/2;
A  = [1 T; 0 1];
B  = [0; T];
nx = 2;   
nu = 1;  
nc = 6;
Q  = eye(nx); 
R  = 0.1*eye(nu);        
[K, Px] = dlqr(A, B, Q, R); 
K       = -K;
Phi     = A + B*K;
N       = 10;
Pc = B'*Px*B + R;
for i  = 1:N-1
    Pc = blkdiag(Pc, B'*Px*B + R);
end

% define constraint_set
distance_error = 6.5; 
acc_bound      = 3;       
speed_bound    = 5;

% joint constraint
%F  = [1/distance_error 0; -1/distance_error 0;0 1/speed_bound; 0 -1/speed_bound; 0 0;0 0]; % Dimension nc x nx
F  = [1/distance_error 0; -1/distance_error 0;0 0; 0 0; 0 0;0 0];
G  = [0; 0; 0; 0; 1/acc_bound; -1/acc_bound];  % Dimension nc x nu

E = eye(nu);
M = [];
for i = 1:N-1
    E = [E zeros(nu)];
    M = blkdiag(M,eye(nu));
end
M     = [zeros((N-1)*nu,nu) M;zeros(nu,N*nu)]; % block-upshift 
Psi   = [Phi B*E;zeros(nu*N, nx) M];
F_bar = [F + G*K G*E]; 

% define true disturbance set and the max. disturbance set
%W true set
theta_EV    = 0;           % rotation angle of \Xi_true_EV (rotation just makes the set more complex)
theta_LV    = 0;           % rotation angle of \Xi_true_LV
Rotation_EV = [cos(theta_EV) -sin(theta_EV); sin(theta_EV) cos(theta_EV)]; % rotation matrix of \Xi_true_EV
Rotation_LV = [cos(theta_LV) -sin(theta_LV); sin(theta_LV) cos(theta_LV)]; % rotation matrix of \Xi_true_LV

Xi_true_EV  = Rotation_EV*Polyhedron([-0.02 -0.15;0.02 -0.15; 0.02 0.15; -0.025 0.15]); % \Xi_true_EV
Xi_true_LV  = Rotation_LV*Polyhedron([-0.02 -0.15;0.02 -0.15; 0.02 0.15; -0.025 0.15]); % \Xi_true_LV
min_u_LV    = -1/15;  % lower bound of acceleration of LV
max_u_LV    = 1/15;   % upper bound of acceleration of LV
U_true_LV   = Polyhedron([1/max_u_LV; 1/min_u_LV], [1; 1]);  % U_true_LV
Wtrue       = Xi_true_EV + (-1*Xi_true_LV) + (-B*U_true_LV);  % true disturbance set of the relative model, W_true
% W set

V_Robust = [-0.2071 -0.5; 0.2071 -0.5; 0.2071 0.5; 0.5 -0.2071; 0.5	0.2071; -0.2071	0.5; -0.5 0.2071; -0.5 -0.2071];
%V_Robust  = [-1 -0.5; 1 -0.5; 1 0.5; -1 0.5];
W = Polyhedron(V_Robust);
W_A = W.A;
W_b = W.b;

W_A    = W_A ./ W_b;
W_b(:) = 1;

W = Polyhedron(W_A, W_b);

% generate sample
Num_off     = 100;
xi_veh0     = Uiniform_Sampling(Xi_true_LV,Num_off);
xi_veh1     = Uiniform_Sampling(Xi_true_EV,Num_off);
u_veh0      = Uiniform_Sampling(U_true_LV,Num_off);
WSample_off = xi_veh1 - xi_veh0 - B*u_veh0;

Num_online     = 100;
xi_veh0        = Uiniform_Sampling(Xi_true_LV,Num_online);
xi_veh1        = Uiniform_Sampling(Xi_true_EV,Num_online);
u_veh0         = Uiniform_Sampling(U_true_LV,Num_online);
WSample_online = xi_veh1 - xi_veh0 - B*u_veh0;

% RPMI Set
% compute MRPI set
epsilon  = 0.15;
Sp       = compute_mrpi_set(Phi, W, epsilon);
SpA      = Sp.A./Sp.b;
Spb      = ones(length(Sp.b), 1);

S       = Polyhedron(SpA, Spb);
S       = minHRep(S); 

h       = Polysupport(S, (F + G*K)');
e_max   = Polysupport(S, (S.A*Phi)');
w_max   = Polysupport(W, (S.A)');

% SAVE to STRUCT
parameters2V = struct;
parameters2V.S     = S;
parameters2V.h     = h;
parameters2V.e_max = e_max;
parameters2V.w_max = w_max;
parameters2V.T = T;
parameters2V.A = A;
parameters2V.B = B;
parameters2V.W = W;
parameters2V.Wtrue = Wtrue;
parameters2V.nx = nx;
parameters2V.nu = nu;
parameters2V.nc = nc;
parameters2V.nv = length(W.b); % W = {w|V*w <= 1}, V: nv x nx
parameters2V.WSample_off = WSample_off;
parameters2V.WSample_online = WSample_online;
parameters2V.num_sample = Num_off;
parameters2V.F = F;
parameters2V.G = G;
parameters2V.K = K;
parameters2V.F_bar = F_bar;
parameters2V.Phi = Phi;
parameters2V.Psi = Psi;
parameters2V.Px = Px;
parameters2V.Pc = Pc;
parameters2V.qa = 0.1;
parameters2V.N = N;
save 'parameters2V.mat' parameters2V;

figure(1)
plot(Wtrue, 'wire', 1, 'edgecolor', 'm', 'linewidth', 2);
hold on
plot(W, 'wire', 1, 'edgecolor', 'r', 'linewidth', 2)

figure(2)
plot(S)


