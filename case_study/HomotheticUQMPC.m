classdef HomotheticUQMPC < handle

    properties (SetAccess = public)
        W;
        Sh;
        WSample_off;
        nx;
        nv;
        num_sample;
        Pa;
        nu;
        N;
        F_bar;
        Psi;
        hh;
        e_max;
        tau_h;
        gamma;
        Px;
        Pc;

        optimizerini;
        optimizerrecursive;
        optimizerWhat;
        optimizerpolysupport;
    end
    
    methods (Access = public)
        function obj = HomotheticUQMPC(parameter)

            obj.W = parameter.W; 
            obj.Sh = parameter.Sh;
            obj.WSample_off = parameter.WSample_off; 
            obj.nx = parameter.nx; 
            obj.nv = parameter.nv; 
            obj.num_sample = parameter.num_sample; 
            obj.Pa = parameter.Pa;
            obj.nu = parameter.nu;
            obj.N = parameter.N;
            obj.F_bar = parameter.F_bar;
            obj.Psi = parameter.Psi;
            obj.hh = parameter.hh;
            obj.e_max = parameter.e_max;
            obj.tau_h = parameter.tau_h;
            obj.gamma = parameter.gamma;
            obj.Px = parameter.Px;
            obj.Pc = parameter.Pc;
        
            obj.optimizerini = obj.optimizer_ini( );
            obj.optimizerrecursive = obj.optimizer_recursive( );
            obj.optimizerWhat = obj.optimizer_W_hat( );
            obj.optimizerpolysupport = obj.optimizer_polysupport( );

        end

        function [theta, rho, v] = initialization_mpc(obj)
            [theta, rho, y] = obj.optimizerini(obj.WSample_off);
            theta = full(theta);
            rho   = full(rho);
            y     = full(y);
            if abs(1-rho) <= 1e-8
                v = obj.WSample_off(:,1);
            else
                v = y/(1-rho);
            end
        end

        function [s, c, alpha] = recursive_mpc(obj, xini, ogmega_max, alpha_max)
            [s, c , alpha] = obj.optimizerrecursive(xini, ogmega_max, alpha_max);
            s      = full(s);
            c      = full(c);
            alpha  = full(alpha);
        end

        function [theta, rho, v, Wb, omega_max, alpha_max] = update_W_hat(obj, theta_old, rho_old, v_old, sample)
            if obj.W.A*sample <= (theta_old + (1 - rho_old)*obj.W.A*v_old)
                theta = theta_old;
                rho   = rho_old;
                v     = v_old;
            else 
                [theta, rho, y] = obj.optimizerWhat(sample, theta_old, rho_old, v_old);
                theta = full(theta);
                rho   = full(rho);
                y     = full(y);
                if abs(1 - rho) <= 1e-8
                    v = sample;
                else
                    v = y/(1 - rho);
                end
            end

            Wb = theta + (1 - rho)*obj.W.A*v;
            WA = (obj.Sh.A)';
            m = size(WA,2);
            omega_max = zeros(m,1);
            for i = 1:m
                item = obj.optimizerpolysupport(Wb, WA(:,i));
                omega_max(i) = WA(:,i)'*full(item);
            end

           alpha_max = max(omega_max)/(1 - obj.gamma);

        end

        function output = optimizer_polysupport(obj)

            opti  = casadi.Opti(); 
            Pb    = opti.parameter(length(obj.W.b), 1);
            n     = size(obj.W.A, 2);
            yalmip('clear')
            z   = opti.variable(n,1);
            y   = opti.parameter(n,1);
            opti.subject_to(obj.W.A*z <= Pb);
            J = y'*z;
            opti.minimize(-J);
            
            options = struct;
            options.ipopt.linear_solver = 'ma57';
            options.ipopt.print_level = 0;
            opti.solver('ipopt', options); 
            output = opti.to_function('f', {Pb, y}, {z});

        end

        function output = optimizer_ini(obj)

            opti   = casadi.Opti(); 
            sample = opti.parameter(obj.nx, obj.num_sample);
            theta  = opti.variable(obj.nv, 1);
            rho    = opti.variable( );
            y      = opti.variable(obj.nx,1);

            for i=1:obj.num_sample
                opti.subject_to(-obj.W.A*y <= theta - obj.W.A*sample(:,i));
            end
            opti.subject_to(obj.W.A*y <= (1 - rho)*obj.W.b);
            opti.subject_to(0 <= rho <=1);
            opti.subject_to(0 <= theta <= rho);
            J = ones(1,obj.nv)*theta + rho;
            opti.minimize(J);

            options = struct;
            options.ipopt.linear_solver = 'ma57';
            options.ipopt.print_level = 0;
            opti.solver('ipopt', options); 
            output = opti.to_function('f', {sample}, {theta, rho, y});
        end

        function output = optimizer_recursive(obj)

            opti  = casadi.Opti( ); 
            s     = opti.variable(obj.nx, 1);
            c     = opti.variable(obj.nu*obj.N, 1);
            alpha = opti.variable(obj.N + 1, 1);
            xini  = opti.parameter(obj.nx, 1);
            omega_max = opti.parameter(length(obj.W.omega_max), 1);
            alpha_max = opti.parameter( );

            opti.subject_to(0 <= alpha);
            opti.subject_to(obj.Sh.A*(xini - s) <= alpha(1)*obj.Sh.b);

            for i = 1:obj.N 
                opti.subject_to(obj.F_bar*obj.Psi^(i-1)*[s;c] <= 1-alpha(i)*obj.hh);
                opti.subject_to(alpha(i)*obj.e_max + omega_max <= alpha(i + 1));
            end
            
            for i = obj.N + 1:obj.tau_h + 1
                opti.subject_to(obj.F_bar*obj.Psi^(i-1)*[s;c] <= 1-(obj.gamma^(i - obj.N-1)*(alpha(obj.N + 1) - alpha_max) + alpha_max)*obj.hh);
            end

            J = s'*obj.Px*s + c'*obj.Pc*c + (alpha - alpha_max)'*obj.Pa*(alpha - alpha_max);
            opti.minimize(J);

            options = struct;
            options.ipopt.linear_solver = 'ma57';
            options.ipopt.print_level = 0;
            opti.solver('ipopt', options); 
            output = opti.to_function('f', {xini, omega_max, alpha_max}, {s, c, alpha});
        end

        function output = optimizer_W_hat(obj)

            opti  = casadi.Opti( ); 
            theta  = opti.variable(obj.nv, 1);
            rho    = opti.variable( );
            y      = opti.variable(obj.nx, 1);

            theta_old  = opti.parameter(obj.nv, 1);
            rho_old    = opti.parameter( );
            v_old      = opti.parameter(obj.nx, 1);
            y_old      = (1 - rho_old)*v_old;
            sample     = opti.parameter(obj.nx, 1);

            opti.subject_to(-obj.W.A*y <= theta - obj.W.A*sample);
            opti.subject_to(-obj.W.A*y <= theta - theta_old - obj.W.A*y_old);
            opti.subject_to(obj.W.A*y <= (1-rho));
            opti.subject_to(0 <= rho <= 1);
            opti.subject_to(0 <= theta <= rho);

            J = ones(1,obj.nv)*theta + rho;

            opti.minimize(J);

            options = struct;
            options.ipopt.linear_solver = 'ma57';
            options.ipopt.print_level = 0;
            opti.solver('ipopt', options); 
            output = opti.to_function('f', {sample, theta_old, rho_old, v_old}, {theta, rho, y});
        end
    end
end