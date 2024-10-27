classdef TradHomotheticMPC < handle

    properties (SetAccess = public)
        W;
        S;
        WSample_off;
        nx;
        nv;
        num_sample;
        nu;
        N;
        F_bar;
        Psi;
        h;
        e_max;
        w_max;
        tau_th;
        Px;
        Pc;
        qa;

        optimizerini;
        optimizerrecursive;
        optimizerWhat;
        optimizerpolysupport;
    end
    
    methods (Access = public)
        function obj = TradHomotheticMPC(parameter)

            obj.W = parameter.W; 
            obj.S = parameter.S;
            obj.WSample_off = parameter.WSample_off; 
            obj.nx = parameter.nx; 
            obj.nv = parameter.nv; 
            obj.num_sample = parameter.num_sample; 
            obj.nu = parameter.nu;
            obj.N = parameter.N;
            obj.F_bar = parameter.F_bar;
            obj.Psi = parameter.Psi;
            obj.h   = parameter.h;
            obj.e_max = parameter.e_max;
            obj.w_max = parameter.w_max;
            obj.tau_th = parameter.tau_th;
            obj.Px = parameter.Px;
            obj.Pc = parameter.Pc;
            obj.qa = parameter.qa;
        
            obj.optimizerrecursive   = obj.optimizer_recursive( );

        end

        function [s, c, alpha] = recursive_mpc(obj, xini)
            [s, c , alpha] = obj.optimizerrecursive(xini);
            s              = full(s);
            c              = full(c);
            alpha          = full(alpha);
        end

        function output = optimizer_recursive(obj)

            opti  = casadi.Opti( ); 
            s     = opti.variable(obj.nx, 1);
            c     = opti.variable(obj.nu*obj.N, 1);
            alpha = opti.variable(obj.N + 1, 1);
            xini  = opti.parameter(obj.nx, 1);

            opti.subject_to(0 <= alpha);
            opti.subject_to(obj.S.A*(xini - s) <= alpha(1)*obj.S.b);

            for i = 1:obj.N 
                opti.subject_to(obj.F_bar*obj.Psi^(i-1)*[s;c] <= 1-alpha(i)*obj.h);
                opti.subject_to(alpha(i)*obj.e_max + obj.w_max <= alpha(i + 1));
            end
            
            for i = obj.N + 1:obj.tau_th + 1
                opti.subject_to(obj.F_bar*obj.Psi^(i-1)*[s;c] <= 1-obj.h);
            end

            J = s'*obj.Px*s + c'*obj.Pc*c + obj.qa*(alpha - 1)'*(alpha - 1);
            opti.minimize(J);

            options = struct;
            options.ipopt.linear_solver = 'ma57';
            options.ipopt.print_level = 0;
            opti.solver('ipopt', options); 
            output = opti.to_function('f', {xini}, {s, c, alpha});
        end

    end
end