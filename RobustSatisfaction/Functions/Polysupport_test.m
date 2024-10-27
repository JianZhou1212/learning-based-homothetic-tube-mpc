function [Y_ut, m, supportvect] = Polysupport_test(P, Y)

    m = size(Y,2);
    Y_ut = Y;
    n = size(P.A, 2);
    yalmip('clear')
    z   = sdpvar(n,1);
    y   = sdpvar(n,1);
    cns = [P.A*z <= P.b];
    obj = y'*z;
    ops = sdpsettings('verbose',0);
    Support_OPT = optimizer(cns, -obj, ops, y, z);
    supportvect = zeros(m,1);

    parfor i = 1:m
        item = Support_OPT(Y(:,i));
        supportvect(i,1) = Y(:,i)'*item;
    end

end

