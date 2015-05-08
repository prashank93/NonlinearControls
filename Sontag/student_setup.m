% Function to setup and perform any one-time computations
% Input parameters
%   x - state of the rocket, x=[y,z,th,psi,dy,dz,dth,dpsi,m]^T
%   consts - structure that contains various system constants
% Output parameters
%   ctrl  -  any student defined control parameters
function ctrl = student_setup(x0, consts)
global done
done = 0;

qr1 = [0.00204073753336390; -0.000908019845082666;94.6797088671666;-0.000565059278701018;-0.000481162874801398;0.00119017653895017;8.89500774975685;0.00379848720593460;-0.000741772844496262;0.529571440663460;1.85372872712185]';
qr = qr1;
qr2 = [0.101641564073661; 96.9396508453469; 11.2229550402918; 0.00107600564745053;-0.00167770724846429;0.00150614659304502;0.000265372697462350;-0.00153529302399248;0.00104849188330871;0.0983966482036966;9.93800621599940]';
%{
if x0(3) > pi/2
    qr = qr2;
end
%}
%{
qr = [0.001 0.001 100 0.001 0.001 0.001 10 0.001 0.001 .5 1];
if x0(3) > pi/2
    qr = [0.0001 0.01 1 0.001 0.0001 0.01 1 0.001 0.001 .1 1];
end
%}
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

ft_0 = (consts.m_nofuel*consts.g)/consts.gamma;

target = [0 consts.L 0 0 0 0 0 0 consts.max.m_fuel];
eq = [0 consts.L  0 0 0 0 0 0 consts.max.m_fuel 1.7168 0];
eq_vars = [y z th ps dy dz dth dpsi m ft T];

A_sym = jacobian(sys,x);
B_sym = jacobian(sys,u);

A = eval(subs(A_sym,eq_vars,eq));
B = eval(subs(B_sym,eq_vars,eq));

Q = diag(abs(qr(1:9)));
R = diag(abs(qr(10:11)));
[K,P] = lqr(A,B,Q,R);
end