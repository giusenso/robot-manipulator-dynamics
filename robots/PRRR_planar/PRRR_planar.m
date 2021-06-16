%%  PRRR planar robot dynamics
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

close all
clear all
clc

%% Symbolic parameters

N   = 4; % number of joints ### ROBOT DATA ###
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
syms I4xx I4xy I4xz I4yx I4yy I4yz I4zx I4zy I4zz
I{1} = [I1xx I1xy I1xz; I1yx I1yy I1yz; I1zx I1zy I1zz];
I{2} = [I2xx I2xy I2xz; I2yx I2yy I2yz; I2zx I2zy I2zz];
I{3} = [I3xx I3xy I3xz; I3yx I3yy I3yz; I3zx I3zy I3zz];
I{4} = [I4xx I4xy I4xz; I4yx I4yy I4yz; I4zx I4zy I4zz];

%% ROBOT DATA
sigma = [1 0 0 0];          % 0 if revolute, 1 if prismatic
syms h;
r{1} = [q(1); h; 0];
r{2} = [l(2)*cos(q(2)); l(2)*sin(q(2)); 0];
r{3} = [l(3)*cos(q(3)); l(3)*sin(q(3)); 0];
r{4} = [l(4)*cos(q(4)); l(4)*sin(q(4)); 0];
gravity = [0; 0; -g0];

%% DYNAMICS COMPUTATION
[rcm0,vcm,w] = moving_frame_gen(q, dq, r, l, d, sigma);
Tcm = kinetic_energy(dq, m, I, vcm, w);
T = sum(Tcm);
M = inertia_matrix(T, q, dq);
[c,christoff,S,X] = centrifugal_coriolis(M, q, dq);
Ucm = potential_energy(m, gravity, rcm0);
U = sum(Ucm);
g = simplify(jacobian(U,q).');

%% Substitutions
% impose equal link length
syms L
[vcm, rcm0, Tcm, M, c, christoff, S, Ucm, g] = ...
    subs_all(l, [L;L;L;L], vcm, rcm0, Tcm, M, c, christoff, S, Ucm, g);


%% Print Dynamics
print_dynamics(sigma, vcm, w, rcm0, Tcm, M, c, christoff, S, Ucm, g)


%