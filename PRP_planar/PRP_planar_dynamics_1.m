%%  PRP planar robot - Dynamics
%   03.2020, @Giuseppe Sensolini
%	Robotics toolbox by Peter Corke required

%clear all
close all
clear all
clc


%% Robot parameters________________________________________________________
disp('::: INPUT DATA :::::::::::::::::::::::::::::::::::');

sigma  = [1 0 1]  % 0 if revolute, 1 if prismatic
N = 3;          % joint number                   
q = sym('q',[N,1]), assume(q,'real');    % joints symbolic variables
l = sym('l',[N,1]), assume(l,'real');    % joints symbolic variables


%% DYNAMICS _______________________________________________________________
disp("DYNAMICS ###########################################################");

syms g0 real;
gravity = [0; 0; -g0];  % gravity vector
z   = [0; 0; 1];
dq  = sym('dq',[N,1]); assume(dq,'real');   % first derivative of joint varaibles
ddq = sym('ddq',[N,1]); assume(ddq,'real'); % second derivative of joint varaibles
dcm = sym('dcm',[N,1]); assume(dcm,'real'); % discance of the CoM from the joint axes
m   = sym('m',[N,1]); assume(m,'real');     % masses of CoM
u   = sym('u',[N,1]); assume(u,'real');     % external forces

% inertia matrices
I = cell(1,N);
syms I1xx I1xy I1xz I1yx I1yy I1yz I1zx I1zy I1zz
syms I2xx I2xy I2xz I2yx I2yy I2yz I2zx I2zy I2zz
syms I3xx I3xy I3xz I3yx I3yy I3yz I3zx I3zy I3zz
I{1} = [I1xx I1xy I1xz; I1yx I1yy I1yz; I1zx I1zy I1zz];
I{2} = [I2xx I2xy I2xz; I2yx I2yy I2yz; I2zx I2zy I2zz];
I{3} = [I3xx I3xy I3xz; I3yx I3yy I3yz; I3zx I3zy I3zz];

% angular and linear velocity
prev_w = zeros(3,1);
prev_v = zeros(3,1);
w = cell(1,N);
v = cell(1,N);
vcm = cell(1,N);
r = cell(1,N);
rcm = cell(1,N);
rcm0 = cell(1,N);

w{1}    = [0; 0; 0];
r{1}    = [0; q(1)+l(1); 0];
rcm{1}  = [0; q(1); 0];
rcm0{1} = rcm{1};
v{1}    = [0; dq(1); 0];
vcm{1}  = [0; dq(1); 0];

w{2}    = dq(2)*z;
r{2}    = [l(2)*cos(q(2)); l(2)*sin(q(2)); 0];
rcm{2}  = [dcm(2)*cos(q(2)); dcm(2)*sin(q(2)); 0];
rcm0{2} = r{1} + rcm{2};
v{2}    = v{1} + cross(w{2},r{2});
vcm{2}  = v{1} + cross(w{2},rcm{2});

w{3}    = w{2};
r{3}    = [-q(3)*sin(q(2)); q(3)*cos(q(2)); 0];
rcm{3}  = [-q(3)*sin(q(2)); q(3)*cos(q(2)); 0];
rcm0{3} = r{1} + r{2} + rcm{3};
v{3}    = v{2} + [-dq(3)*sin(q(2)); dq(3)*cos(q(2)); 0] + cross(w{3},r{3});
vcm{3}  = v{2} + [-dq(3)*sin(q(2)); dq(3)*cos(q(2)); 0] + cross(w{3},rcm{3});

disp("::: KINETIC ENERGY ::::::::::::::::::::::::::::::::::::::::::::::::::");
[T, Tcm] = kinetic_energy(q, dq, m, I, vcm, w)

disp("::: INERTIA MATRIX ::::::::::::::::::::::::::::::::::::::::::::::::::");
M = inertia_matrix(T, q, dq)
dM = inertia_diff(M, q, dq);

disp("::: CORIOLIS AND CENTRIFUGAL ::::::::::::::::::::::::::::::::::::::::");
[c, christoff, S, X] = centrifugal_coriolis(M, q, dq)

disp("::: GRAVITY :::::::::::::::::::::::::::::::::::::::::::::::::::::::::");
[U, Ucm] = potential_energy(q, m, dcm, gravity, rcm0)
g = simplify(jacobian(U,q).')


%% Parametrization

params = [  m(1) + m(2) + m(3);
            I{2}(3,3) + I{3}(3,3) + m(2)*(dcm(2)^2);
            m(3);
            m(2)*dcm(2)]

[M_a, a] = parametrization(M, params);
[c_a, a] = parametrization(c, params);
[g_a, a] = parametrization(g, params);

Y = simplify(jacobian(M_a*ddq+c_a,a))

















%