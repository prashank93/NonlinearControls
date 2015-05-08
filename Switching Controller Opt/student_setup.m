% Function to setup and perform any one-time computations
% Input parameters
%   x - state of the rocket, x=[y,z,th,psi,dy,dz,dth,dpsi,m]^T
%   consts - structure that contains various system constants
% Output parameters
%   ctrl  -  any student defined control parameters
function ctrl = student_setup(x0, consts, qr1, qr2)
global done
done = 0;

% qr = [0.001 0.001 100 0.001 0.001 0.001 10 0.001 0.001 .5 2];
qr = qr1;
if x0(3) > 0.75*pi && x0(2) < 400
    %     qr = [0.1 100 10 0.000001 0.001 0.001 .001 0.001 0.001 .1 10];
    qr = qr2;
end
[~,P] = lqrmaker(qr);
K = [ -0.7544 1.1544 0 0 -5.5728 8.5728 0 0 -0.1575;
    0.1306    0  -37.6338   72.8933    1.6009    0  -48.6607   30.3236    0];
% [K] = lqrmaker(qr2);
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 real

f_vec = [   x5
    x6
    x7
    x8
    0
    -consts.g
    0
    0
    0] ;

g_vec = [0,    0 ;
    0,    0 ;
    0,    0 ;
    0,    0 ;
    -consts.gamma*sin(x4+x3)/x9,    0 ;
    consts.gamma*cos(x4+x3)/x9,    0 ;
    -consts.L*consts.gamma*sin(x4)/consts.J,    0 ;
    0, 1/consts.JT ;
    -1,    0] ;

x = [x1 x2 x3 x4 x5 x6 x7 x8 x9]';
V = x'*P*x;
LfV(x1,x2,x3,x4,x5,x6,x7,x8,x9) = jacobian(V,x)*f_vec;
LgV(x1,x2,x3,x4,x5,x6,x7,x8,x9) = jacobian(V,x)*g_vec;
ctrl.LfVfun = matlabFunction(LfV);
ctrl.LgVfun = matlabFunction(LgV);
ctrl.K = K;
end

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

ft_0 = ((consts.m_nofuel+consts.max.m_fuel)*consts.g)/consts.gamma;

target = [0 consts.L 0 0 0 0 0 0 consts.max.m_fuel];
%1.4715
eq = [0 consts.L 0 0 0 0 0 0 consts.max.m_fuel 1.7168 0];
eq_vars = [y z th ps dy dz dth dpsi m ft T];

A_sym = jacobian(sys,x);
B_sym = jacobian(sys,u);

A = eval(subs(A_sym,eq_vars,eq));
B = eval(subs(B_sym,eq_vars,eq));

Q = diag(abs(qr(1:9)));
R = diag(abs(qr(10:11)));
[K,P] = lqr(A,B,Q,R);
end