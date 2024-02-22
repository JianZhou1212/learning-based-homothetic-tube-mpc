function [P,p] = findMRPIset(Phi,C,c,Pw,pw,maxiters,tol,lp_params,print_flag)

if nargin < 6 || isempty(maxiters)
    maxiters = 20;
end
if nargin < 7 || isempty(tol)
    tol = 1e-7;
end
if nargin < 8 || isempty(lp_params)
    lp_params.OutputFlag = 0; % Gurobi runs silently
end
if nargin < 9 || isempty(print_flag)
    print_flag = 0; % Run silently
end

[nC,nx] = size(C);
[nW,nw] = size(Pw);

C_k = C;
e_k = zeros(nC,1);
P = C;
p = c;

lp_e.rhs = pw;
lp_e.A = sparse(Pw);

stop_flag = 0;
infeas_flag = 0;
k = 1;
while ~stop_flag && ~infeas_flag && k < maxiters

    C_knext = C_k*Phi;
    e_knext = zeros(nC,1);

    stop_flag = 1;

    yalmip('clear')
    z=sdpvar(nx,1);
    y=sdpvar(nx,1);
    cns=[];
    cns=[cns,lp_e.A*z<=lp_e.rhs];
    obj=y'*z;
    ops = sdpsettings('verbose',0);
    OPTe = optimizer(cns,-obj, ops, y, z);

    for i = 1:nC
        item=OPTe(C_k(i,:)');
        e_knext(i) = C_k(i,:)*item + e_k(i);

        if stop_flag
            lp_t.rhs = p;
            lp_t.A = sparse(P);
            yalmip('clear')
            z=sdpvar(nx,1);
            cns=[];
            cns=[cns,lp_t.A*z<=lp_t.rhs];
            objt=C_knext(i,:)*z;
            ops = sdpsettings('verbose',0);
            OPTt=optimize(cns,-objt,ops);

            if ~OPTt.problem
                infeas_flag = 1;
                stop_flag = 0;
            else
                t = value(objt) + e_knext(i) - c(i);
                if print_flag > 1
                    fprintf(1,'%d %d %.4e\n', k, i, t);
                end
                if t > 0
                    stop_flag = 0;
                end
            end
        end
    end
  
    if ~stop_flag
        P = [P; C_knext]; %#ok<AGROW>
        p = [p; c-e_knext]; %#ok<AGROW>
        C_k = C_knext;
        e_k = e_knext;
        k = k+1;
    end
end
end
