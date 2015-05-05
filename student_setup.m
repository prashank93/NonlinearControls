% Function to setup and perform any one-time computations
% Input parameters
%   x - state of the rocket, x=[y,z,th,psi,dy,dz,dth,dpsi,m]^T
%   consts - structure that contains various system constants
% Output parameters
%   ctrl  -  any student defined control parameters
function ctrl = student_setup(x0, consts)
% qr = [100 1 100 1 50 1 100 1 1 1 1];
qr = ones(1,11);
[~,P] = lqrmaker(qr);

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
% ctrl.K = K;
end