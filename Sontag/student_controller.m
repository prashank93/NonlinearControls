% Student supplied function to compute the control input at each instant of time
% Input parameters:
%   t  -  Time (in seconds) since start of simulation
%   x - state of the rocket, x=[y,z,th,psi,dy,dz,dth,dpsi,m]^T
%   consts - structure that contains various system constants
%   ctrl  -  any student defined control parameters
% Output parameters:
%   u  -  [thrust; torque] - two inputs to the rocket
function u = student_controller(t, x, consts, ctrl)
global done

% LfV = eval(subs(ctrl.LfV,sym_vec,x));
% LgV = eval(subs(ctrl.LgV,sym_vec,x));


if done == 0 && x(2) > 200
    x(2) = x(2) - 200;
end


LfV = ctrl.LfVfun(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9));
LgV = ctrl.LgVfun(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9));
% a = LfV;
% b = LgV*LgV';
% if LgV ~= 0
%     u = -((a + sqrt(a^2 + b^4))/b)*LgV';
% else
%     u = [0;0];
% end

if all(LgV ~= 0)
   %{ 
    u(1,1) = -((LfV + sqrt(LfV^2 + LgV(1)^4))/LgV(1));
    u(2,1) = -((LfV + sqrt(LfV^2 + LgV(2)^4))/LgV(2));
%}
    u = -LgV'*(LfV + sqrt(LfV^2 + (LgV*LgV')^2))/(LgV*LgV');
else
    u = [0;0];
end
% && abs(x(5)) < 5
% if abs(x(3)) < 0.05 && abs(x(7)) < .1 && abs(x(2)) < 100 || done == 1
if abs(x(3)) < 0.05 && abs(x(7)) < .1 && abs(x(2)) < 40 || done == 1
    u = -ctrl.K*x;
    done = 1;
end
end