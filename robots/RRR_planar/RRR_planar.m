%%  RRR planar robot dynamics
%   Developed by Giuseppe Sensolini (github: @giusenso)

%{
This script can be adapted to any robot manipulator.
Variables to be adjusted:
    - N (number of joints)
    - sigma vector (0=revolute joint, 1=prismatic joint)
    - r{i} (vector going from the from joint i to i+1 in the world reference frame)
    - gravity vector
    - eventually add other inertial symbolic terms if needed
%}

clear all
close all
clc

%% Symbolic parameters

N   = 3; % number of joints ### ROBOT DATA ###
q   = sym('q',[N,1]); assume(q,'real');    	% joints symbolic variables
l   = sym('l',[N,1]); assume(l,'real');   	% joints symbolic variables
d   = sym('d',[N,1]); assume(d,'real'); % distance of the CoM from the joint axes
dq  = sym('dq',[N,1]); assume(dq,'real');   % first derivative of joint varaibles
ddq = sym('ddq',[N,1]); assume(ddq,'real'); % second derivative of joint varaibles
m   = sym('m',[N,1]); assume(m,'real');     % masses of CoM
u   = sym('u',[N,1]); assume(u,'real');     % external forces
r = cell(1,N);  % vector going from frame i-1 to frame i expressed in frame 0
syms g0 real;

% inertia matrices
I = cell(1,N);
syms I1xx I1xy I1xz I1yx I1yy I1yz I1zx I1zy I1zz
syms I2xx I2xy I2xz I2yx I2yy I2yz I2zx I2zy I2zz
syms I3xx I3xy I3xz I3yx I3yy I3yz I3zx I3zy I3zz
I{1} = [I1xx I1xy I1xz; I1yx I1yy I1yz; I1zx I1zy I1zz];
I{2} = [I2xx I2xy I2xz; I2yx I2yy I2yz; I2zx I2zy I2zz];
I{3} = [I3xx I3xy I3xz; I3yx I3yy I3yz; I3zx I3zy I3zz];

%% ROBOT DATA
sigma = [0 0 0];        % 0 if revolute, 1 if prismatic
r{1} = [l(1)*cos(q(1)); l(1)*sin(q(1)); 0];
r{2} = [l(2)*cos(q(1)+q(2)); l(1)*sin(q(1)+q(2)); 0];
r{3} = [l(3)*cos(q(1)+q(2)+q(3)); l(3)*sin(q(1)+q(2)+q(3)); 0];
gravity = [0; -g0; 0];

%% DYNAMICS COMPUTATION
[rcm0,vcm,w] = moving_frame(q, dq, r, l, d, sigma);
Tcm = kinetic_energy(dq, m, I, vcm, w);
T = simplify(sum(Tcm));
M = inertia_matrix(T, q, dq);
[c,christoff,S,X] = centrifugal_coriolis(M, q, dq);
Ucm = potential_energy(m, gravity, rcm0);
U = simplify(sum(Ucm));
g = simplify(jacobian(U,q).');

%% Print Dynamics
print_dynamics(sigma, vcm, w, rcm0, Tcm, M, c, christoff, S, Ucm, g)


%{
%% kinemetics
f = r{1} + r{2} + r{3}
J = subs(jacobian(f(1:2),q), {l(1),l(2),l(3), d(1),d(2),d(3), m(1),m(2),m(3)}, {.5,.4,.3, .2,.1,.2, 5,3,2})
U = subs(U, {g0, l(1),l(2),l(3), d(1),d(2),d(3), m(1),m(2),m(3)}, {9.81, .5,.4,.3, .2,.1,.2, 5,3,2})
g_ = double(subs(jacobian(U,q).', {q(1),q(2),q(3)}, {pi,0,-pi/2}))

%% Reduced gradient

% cartesian trajectory planning
q0 = [pi; 0; -pi/2]     % joint configuration
vd = [0; 0]                 % desired e-e velocity
J0 = double(subs(J, q, q0))         

%%
[Ja, Jb, x] = matrix_decoupling(J0)
[qadot, qbdot] = reduced_gradient(Ja, Jb, U, vd, q);
qdot = matrix_coupling(qadot, qbdot, x);
qdot = double(subs(qdot, q, q0))
%}





%