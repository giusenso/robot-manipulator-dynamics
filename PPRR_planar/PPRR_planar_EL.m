%%  PPRR planar robot - Dynamics
%   04.2020, @Giuseppe Sensolini
%	Robotics toolbox by Peter Corke required

%clear all
close all
clear all
clc


%% Robot parameters________________________________________________________
disp('::: INPUT DATA :::::::::::::::::::::::::::::::::::');

syms a alpha d theta
N = 4;      % joint number 
sigma  = [1 1 0 0];  % 0 if revolute, 1 if prismatic
q = sym('q',[N,1]), assume(q,'real');    % joints symbolic variables
l = sym('l',[N,1]), assume(l,'real');    % joints symbolic variables


%% DYNAMICS _______________________________________________________________
disp("DYNAMICS ###########################################################");

syms g0 real;
gravity = [0; -g0; 0];  % gravity vector
z   = [0; 0; 1];
dq  = sym('dq',[N,1]); assume(dq,'real');   % first derivative of joint varaibles
ddq = sym('ddq',[N,1]); assume(ddq,'real'); % second derivative of joint varaibles
d = sym('d',[N,1]); assume(d,'real'); % discance of the CoM from the joint axes
m   = sym('m',[N,1]); assume(m,'real');     % masses of CoM
u   = sym('u',[N,1]); assume(u,'real');     % external forces

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

% angular and linear velocity
r = cell(1,N);
rcm = cell(1,N);
rcm0 = cell(1,N);
w = cell(1,N);
v = cell(1,N);
vcm = cell(1,N);


r{1} = [q(1);0;0];
rcm{1} = [q(1)-d(1);0;0];
rcm0{1} = rcm{1};
w{1} = [0;0;0];
v{1} = [dq(1);0;0];
vcm{1} = [dq(1);0;0];

r{2} = [0; q(2); 0];
rcm{2} = [0; q(2)-d(2); 0];
rcm0{2} = r{1} + rcm{2};
w{2} = [0; 0; 0];
v{2} = v{1} + [0; dq(2); 0];
vcm{2} = v{1} + [0; dq(2); 0];

r{3} = [l(3)*cos(q(3)); l(3)*sin(q(3)); 0];
rcm{3} = [d(3)*cos(q(3)); d(3)*sin(q(3)); 0];
rcm0{3} = r{1}+r{2} + rcm{3};
w{3} = [0; 0; dq(3)];
v{3} = v{2} + cross(w{3},r{3});
vcm{3} = v{2} + cross(w{3},rcm{3});

r{4} = [l(4)*cos(q(3)+q(4)); l(4)*sin(q(3)+q(4)); 0];
rcm{4} = [d(4)*cos(q(3)+q(4)); d(4)*sin(q(3)+q(4)); 0];
rcm0{4} = r{1}+r{2}+r{3} + rcm{4};
w{4} = w{3}+[0; 0; dq(4)];
v{4} = v{3} + cross(w{4},r{4});
vcm{4} = v{3} + cross(w{4},rcm{4});

disp("::: KINETIC ENERGY ::::::::::::::::::::::::::::::::::::::::::::::::::");
[T, Tcm] = kinetic_energy(q, dq, m, I, vcm, w)

disp("::: INERTIA MATRIX ::::::::::::::::::::::::::::::::::::::::::::::::::");
M = inertia_matrix(q, dq, T)

disp("::: CORIOLIS AND CENTRIFUGAL ::::::::::::::::::::::::::::::::::::::::");
[c, christoff] = centrifugal_coriolis(q, dq, M)

disp("::: GRAVITY :::::::::::::::::::::::::::::::::::::::::::::::::::::::::");
[U, Ucm] = potential_energy(q, m, d, gravity, rcm0)
g = simplify(jacobian(U,q).')


%% Parametrization
params = [  m(1) + m(2) + m(3) + m(4);
            m(2) + m(3) + m(4);
            I{3}(3,3) + m(3)*d(3)^2 + I{4}(3,3) + m(4)*d(4)^2 + m(4)*l(3)^2;
            I{4}(3,3) + m(4)*d(4)^2;
            m(4)*d(4);
            m(3)*d(3) + m(4)*l(3)]

[M_a, a] = parametrization(M, params);
[c_a, a] = parametrization(c, params);
[g_a, a] = parametrization(g, params);

Y = simplify( jacobian(M_a*ddq+c_a, a) )






















%