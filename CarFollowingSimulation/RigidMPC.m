classdef RigidMPC < handle
    
    properties (SetAccess = public)

        W;
        S;
        WSample_off;
        N;
        nx;
        nu;
        nc;
        num_sample;
        F;
        G;
        K;
        Px;
        Pc;
        F_bar;
        Psi;
        Phi;
        num_half_space_S;
        tau_k;

        optimizerini;
        optimizerrecursive;
        optimizerWhat;
    end
    
    methods (Access = public)
        function obj = RigidMPC(parameters)

            obj.W = parameters.W;
            obj.S = parameters.S;
            obj.WSample_off = parameters.WSample_off; 
            obj.N = parameters.N;
            obj.nx = parameters.nx;
            obj.nu = parameters.nu;
            obj.nc = parameters.nc;
            obj.num_sample = parameters.num_sample;
            obj.F = parameters.F;
            obj.G = parameters.G;
            obj.K = parameters.K;
            obj.Px = parameters.Px;
            obj.Pc = parameters.Pc;
            obj.F_bar = parameters.F_bar;
            obj.Psi = parameters.Psi;
            obj.Phi = parameters.Phi;
            obj.num_half_space_S = parameters.num_half_space_S;
            obj.tau_k = parameters.tau_k;

            obj.optimizerini = obj.optimizer_ini( );
            obj.optimizerrecursive = obj.optimizer_recursive( );
            obj.optimizerWhat = obj.optimizer_W_hat( );
        end

        function [rho, v] = initialization_mpc(obj)
            [beta, y] = obj.optimizerini(obj.WSample_off);
            beta = full(beta);
            y    = full(y);
            rho  = 1 - beta;
            v    = y/beta;
        end

        function [s, c] = recursive_mpc(obj, xini, Sk, hk)

            S_A = Sk.A;
            S_b = Sk.b;
            
            [s, c] = obj.optimizerrecursive(xini, hk, S_A, S_b);
            s = full(s);
            c = full(c);
        end

        function [rho, v] = update_W_hat(obj, rho_old, v_old, sample)
            if obj.W.A*sample <= (rho_old)*obj.W.b + (1 - rho_old)*obj.W.A*v_old
                rho = rho_old;
                v   = v_old;
            else
                [beta, y] = obj.optimizerWhat(sample, rho_old, v_old);
                beta = full(beta);
                y    = full(y);

                rho  = 1 - beta;
                v      = y/beta;
            end
        end
            
        
        function output = optimizer_ini(obj)
            opti   = casadi.Opti(); 
            sample = opti.parameter(obj.nx, obj.num_sample);
            beta   = opti.variable( );
            y      = opti.variable(obj.nx, 1);

            for i = 1:obj.num_sample
                opti.subject_to(-obj.W.A*y <= (1 - beta) - obj.W.A*sample(:,i));
            end
            opti.subject_to(obj.W.A*y <= beta);
            opti.subject_to(0 <= beta <=1);

            J = -beta;
            opti.minimize(J);

            options = struct;
            options.ipopt.linear_solver = 'ma57';
            options.ipopt.print_level = 0;
            opti.solver('ipopt', options); 
            output = opti.to_function('f', {sample}, {beta, y});
        end
        
        function output = optimizer_recursive(obj)
            opti  = casadi.Opti( );
            s     = opti.variable(obj.nx, 1);
            c     = opti.variable(obj.nu*obj.N, 1);
            xini  = opti.parameter(obj.nx, 1);
            hk    = opti.parameter(obj.nc, 1);
            S_A   = opti.parameter(obj.num_half_space_S, obj.nx);
            S_b   = opti.parameter(obj.num_half_space_S, 1);
            J     = s'*obj.Px*s + c'*obj.Pc*c;
            
            opti.minimize(J);
            opti.subject_to(S_A*(xini - s) <= S_b);
            for i = 0:1:obj.tau_k
                opti.subject_to(obj.F_bar*(obj.Psi^i)*[s; c] <= ones(obj.nc, 1) - hk);
            end

            options = struct;
            options.ipopt.linear_solver = 'ma57';
            options.ipopt.max_iter      = 10000;
            options.ipopt.print_level   = 0;
            opti.solver('ipopt', options); 
            output = opti.to_function('f', {xini, hk, S_A, S_b}, {s, c});
        end

        function output = optimizer_W_hat(obj)
            opti = casadi.Opti( );
            beta = opti.variable( );
            y    = opti.variable(obj.nx, 1);

            rho_old = opti.parameter( );
            v_old   = opti.parameter(obj.nx, 1);
            y_old   = (1 - rho_old)*v_old;
            sample  = opti.parameter(obj.nx, 1);

            opti.minimize(-beta);
            opti.subject_to(0 <= beta <= 1);
            opti.subject_to(obj.W.A*y - beta*obj.W.b <= 0);
            opti.subject_to(rho_old*obj.W.b + obj.W.A*y_old <= (1 - beta)*obj.W.b + obj.W.A*y);
            opti.subject_to(obj.W.A*sample <= (1 - beta)*obj.W.b + obj.W.A*y);

            options = struct;
            options.ipopt.linear_solver = 'ma57';
            options.ipopt.print_level = 0;
            opti.solver('ipopt', options); 
            output = opti.to_function('f', {sample, rho_old, v_old}, {beta, y});
        end


    end
end
