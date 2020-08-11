%%  PRP planar robot - Dynamics
%   04.2020, @Giuseppe Sensolini

clear all
close all
clc

disp("### DYNAMICS ########################################################");

%% Robot parameters________________________________________________________

N = 3  % joint number
                 
q   = sym('q',[N,1]); assume(q,'real');    	% joints symbolic variables
l   = sym('l',[N,1]); assume(l,'real');   	% joints symbolic variables
d   = sym('d',[N,1]); assume(d,'real');   	% distance of the CoM from the joint axes
dq  = sym('dq',[N,1]); assume(dq,'real');   % first derivative of joint varaibles
ddq = sym('ddq',[N,1]); assume(ddq,'real'); % second derivative of joint varaibles
m   = sym('m',[N,1]); assume(m,'real');     % masses of CoM
u   = sym('u',[N,1]); assume(u,'real');     % external forces
r   = cell(1,N);    % vector going from frame i-1 to frame i
syms g0

% inertia matrices
I = cell(1,N);
syms I1xx I1yy I1zz  I2xx I2yy I2zz  I3xx I3yy I3zz  I4xx I4yy I4zz
I{1} = [I1xx 0 0; 0 I1yy 0; 0 0 I1zz];
I{2} = [I2xx 0 0; 0 I2yy 0; 0 0 I2zz];
I{3} = [I3xx 0 0; 0 I3yy 0; 0 0 I3zz];
I{4} = [I4xx 0 0; 0 I4yy 0; 0 0 I4zz];


%### USER DATA HERE! ###
gravity = [0; 0; -g0];  % gravity vector
sigma = [1 1 0 0]       % 0 if revolute, 1 if prismatic 
r{1} = [0; q(1); 0];
r{2} = [l(2)*cos(q(2)); l(2)*sin(q(2)); 0];
r{3} = [-q(3)*sin(q(2)); q(3)*cos(q(2)); 0];
% #######################

disp("::: MOVING FRAME ALGO :::::::::::::::::::::::::::::::::::::::::::::::");
[vcm, w, rcm0] = moving_frame_algo(q, dq, r, l, d, sigma);

disp("::: KINETIC ENERGY ::::::::::::::::::::::::::::::::::::::::::::::::::");
[T, Tcm] = kinetic_energy(q, dq, m, I, vcm, w)

disp("::: INERTIA MATRIX ::::::::::::::::::::::::::::::::::::::::::::::::::");
M = inertia_matrix(T, q, dq)
dM = inertia_diff(M, q, dq);

disp("::: CORIOLIS AND CENTRIFUGAL ::::::::::::::::::::::::::::::::::::::::");
[c, christoff, S, X] = centrifugal_coriolis(M, q, dq)

disp("::: GRAVITY :::::::::::::::::::::::::::::::::::::::::::::::::::::::::");
[U, Ucm] = potential_energy(q, m, d, gravity, rcm0)
g = simplify(jacobian(U,q).')




