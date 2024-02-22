clc
clear 
close all
addpath("Functions");

% system parameters
T     = 1/2;
Vec_A = [1 T; 0 1];
Vec_B = [T^2/2; T];
nx    = 4;  
nu    = 2;   
nc    = 8;
A     = blkdiag(Vec_A, Vec_A);
B     = [Vec_B zeros(2, 1); -Vec_B Vec_B];
Q     = eye(nx); 
R     = 0.1*eye(nu);       
[K, Px] = dlqr(A, B, Q, R);
K       = -K;
Phi     = A + B*K;
N       = 10;
Pc      = B'*Px*B + R;
for i   = 1:N-1
    Pc  = blkdiag(Pc, B'*Px*B + R);
end

% define constraint_set
distance_error = 15; 
acc_bound      = 2.5;       
speed_bound    = 5;

% NOTE: Fx and fx are defined for computing homothetic tube
Fx             = 9*[4/distance_error 0 0 0;
                    -4/distance_error 0 0 0;
                    0 1/(speed_bound) 0 0;
                    0 -1/(speed_bound) 0 0;
                    0 0 4/distance_error 0;
                    0 0 -4/distance_error 0;
                    0  0 0 1/(speed_bound);
                    0  0 0 -1/(speed_bound)]; 
fx = [1;1;1;1;1;1;1;1];
F  = [1/distance_error 0 0 0;
      -1/distance_error 0 0 0;
      0 0 1/distance_error 0;
      0 0 -1/distance_error 0;
      0 0 0 0;
      0 0 0 0;
      0 0 0 0;
      0 0 0 0]; 
G = [0 0; 
     0 0;
     0 0;
     0 0;
     1/acc_bound 0;
     -1/acc_bound 0;
     0 1/acc_bound;
     0 -1/acc_bound];  

E = eye(nu);
M = [ ];
for i = 1:N-1
    E = [E zeros(nu)];
    M = blkdiag(M,eye(nu));
end
M     = [zeros((N-1)*nu, nu) M;zeros(nu, N*nu)];
Psi   = [Phi B*E;zeros(nu*N, nx) M];
F_bar = [F + G*K G*E]; 

% define true disturbance set and the max. disturbance set
%W true set
theta_EV    = pi/40;           % rotation angle of \Xi_true_EV (rotation just makes the set more complex)
theta_LV    = pi/40;           % rotation angle of \Xi_true_LV
Rotation_EV = [cos(theta_EV) -sin(theta_EV); sin(theta_EV) cos(theta_EV)]; % rotation matrix of \Xi_true_EV
Rotation_LV = [cos(theta_LV) -sin(theta_LV); sin(theta_LV) cos(theta_LV)]; % rotation matrix of \Xi_true_LV
Xi_true_EV  = Rotation_EV*Polyhedron([-0.06 -0.015;0.06 -0.015; 0.01 0.025; -0.01 0.025]); % \Xi_true_EV
Xi_true_LV  = Rotation_LV*Polyhedron([-0.06 0.015;0.06 0.015; 0.01 -0.025; -0.01 -0.025]); % \Xi_true_LV
min_u_LV    = -1/20;  % lower bound of acceleration of LV
max_u_LV    = 1/16;   % upper bound of acceleration of LV
U_true_LV   = Polyhedron([1/max_u_LV; 1/min_u_LV], [1; 1]);       % U_true_LV
WtrueMPT1   = Xi_true_EV + (-1*Xi_true_LV) + (-Vec_B*U_true_LV);  % true disturbance set of the relative model, W_true
WtrueMPT2   = Xi_true_EV + (-1*Xi_true_EV);                       % true disturbance set of the relative model, W_true
Wtrue.A     = blkdiag(WtrueMPT1.A, WtrueMPT2.A);
Wtrue.b     = [WtrueMPT1.b; WtrueMPT2.b];

% W set
[W.A, theta,imax] = generate_polytope(nx,1);
W.A = W.A*2.5;
W.b = ones(size(W.A,1),1);
%W   = Polyhedron(W.A, W.b);

% generate sample
Num_off = 100;
xi_veh0 = Uiniform_Sampling(Xi_true_LV, Num_off);
xi_veh1 = Uiniform_Sampling(Xi_true_EV, Num_off);
xi_veh2 = Uiniform_Sampling(Xi_true_EV, Num_off);
u_veh0  = Uiniform_Sampling(U_true_LV,  Num_off);
WSample_off = [xi_veh1 - xi_veh0 - Vec_B*u_veh0
               xi_veh2 - xi_veh1];

Num_online = 100;
xi_veh0    = Uiniform_Sampling(Xi_true_LV, Num_online);
xi_veh1    = Uiniform_Sampling(Xi_true_EV, Num_online);
xi_veh2    = Uiniform_Sampling(Xi_true_EV, Num_online);
u_veh0     = Uiniform_Sampling(U_true_LV,  Num_online);
WSample_online = [xi_veh1 - xi_veh0 - Vec_B*u_veh0
                  xi_veh2 - xi_veh1];

% SAVE to STRUCT
parameters3V = struct;
parameters3V.T = T;
parameters3V.A = A;
parameters3V.B = B;
parameters3V.W = W;
parameters3V.Wtrue = Wtrue;
parameters3V.nx = nx;
parameters3V.nu = nu;
parameters3V.nc = nc;
parameters3V.nv = length(W.b); % W = {w|V*w <= 1}, V: nv x nx
parameters3V.WSample_off = WSample_off;
parameters3V.WSample_online = WSample_online;
parameters3V.num_sample = Num_off;
parameters3V.Fx = Fx;
parameters3V.fx = fx;
parameters3V.F = F;
parameters3V.G = G;
parameters3V.K = K;
parameters3V.F_bar = F_bar;
parameters3V.Phi = Phi;
parameters3V.Psi = Psi;
parameters3V.Px = Px;
parameters3V.Pc = Pc;
parameters3V.N = N;
save 'parameters3V.mat' parameters3V;
