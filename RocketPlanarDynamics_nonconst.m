clear

% symbolic quantities
syms y z th psi real
syms dy dz dth dpsi real
syms ddy ddz ddth ddpsi real

% Inputs
syms f_T tau

% symbolic constants
syms m J JT g L gamma real

% Configuration coordinates and velocities
q = [y;z;th;psi]; 
dq = [dy;dz;dth;dpsi];
ddq = [ddy;ddz;ddth;ddpsi];
u = [f_T ; tau] ;

% Kinetic Energy
KE = 1/2*m*(dy^2+dz^2) + 1/2*J*dth^2 + 1/2*JT*dpsi^2 ;

% Potential Energy
PE = m*g*z ;

% Lagrangian
Lag = KE - PE ;

% Lagrangian equations of motion
LHS = jacobian(jacobian(Lag, dq)', [q;dq])*[dq;ddq]  -  jacobian(Lag, q)' ;
D = jacobian(LHS, ddq) ;
H = simple(LHS - D*ddq) ;

% Input matrix
p_Thruster = [y;z] + L*[sin(th); -cos(th)] ; % Thruster position
B_Thruster = jacobian(p_Thruster,q)' * [-sin(th+psi) ; cos(th+psi)]*gamma ; % Force at thruster position is [-f_T*sin(psi) ; f_T*cos(psi)]
B_vectoring = jacobian(psi, q)' ;
B = simple([B_Thruster  B_vectoring]) ;

% Nonlinear Model (Mass non-constant)
syms dm ddm real ;
x = [q; dq; m] ;

% Drift and Control vector fields
% Note that last-line of f_vec, g_vec is to add dm/dt = -f_T
f_vec = [dq ; 
         -inv(D)*H ;
         0] ;
g_vec = [zeros(length(q), 2) ;
         inv(D)*B ;
         -1 0] ;

 
% Linearization about hover:
consts = get_consts() ;
m_const = consts.m_nofuel+consts.max.m_fuel ;
% linearize about x_lin, u_lin (corresponding to hover)
x_lin = [0;consts.L;0;0 ; 0;0;0;0; m_const] ;
u_lin = [m_const*g;0] ;
A_linear = subs( subs(jacobian(f_vec+g_vec*u, x), x, x_lin), u, u_lin) ;
B_linear = subs(jacobian(f_vec+g_vec*u, u), x, x_lin) ;

A = double(subs(A_linear, {m,J,JT,g,L,gamma}, {consts.m_nofuel+consts.max.m_fuel,consts.J,consts.JT,consts.g,consts.L,consts.gamma})) ;
B = double(subs(B_linear, {m,J,JT,g,L,gamma}, {consts.m_nofuel+consts.max.m_fuel,consts.J,consts.JT,consts.g,consts.L,consts.gamma})) ;

% Check controllability
size(A)
rank(ctrb(A,B))