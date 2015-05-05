function [K,P] = lqrmaker(qr)
consts = get_consts();
global target
syms y z th ps dy dz dth dpsi m ft T real

f_vec = [   dy
            dz
            dth
            dpsi
            0
            -consts.g
            0
            0
            0] ;

g_vec = [0,    0 ;
        0,    0 ;
        0,    0 ;
        0,    0 ;
        -consts.gamma*sin(ps+th)/m,    0 ;
        consts.gamma*cos(ps+th)/m,    0 ;
        -consts.L*consts.gamma*sin(ps)/consts.J,    0 ;
        0, 1/consts.JT ;
        -1,    0] ;

x = [y z th ps dy dz dth dpsi m];
u = [ft; T];

sys = f_vec + g_vec*u;

% consts.m_nofuel
% consts.max.m_fuel
% 1.4715

ft_0 = (consts.m_nofuel*consts.g)/consts.gamma;

target = [0 consts.L 0 0 0 0 0 0 consts.max.m_fuel];
eq = [0 consts.L 0 0 0 0 0 0 consts.max.m_fuel 1.4715 0];
eq_vars = [y z th ps dy dz dth dpsi m ft T];

A_sym = jacobian(sys,x);
B_sym = jacobian(sys,u);

A = eval(subs(A_sym,eq_vars,eq));
B = eval(subs(B_sym,eq_vars,eq));

Control = ctrb(A,B);
% fprintf('Rank = %g\n',rank(Control));
% 
% Q = diag([1 1 1 1 10 10 10 10 1]);
% R = diag([.01 .1]);
Q = diag(abs(qr(1:9)));
R = diag(abs(qr(10:11)));
[K,P] = lqr(A,B,Q,R);

% P2 = lyap(A,Q)
% K = 0;
end

