%%  PRR planar robot - Dynamics
%   04.2020, @Giuseppe Sensolini

%clear all
close all
clear all
clc


%% Robot parameters________________________________________________________
disp('::: INPUT DATA :::::::::::::::::::::::::::::::::::');

syms a alpha d theta
sigma  = [1 0 0]  % 0 if revolute, 1 if prismatic
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
syms I1xx I1yy I1zz
syms I2xx I2yy I2zz
syms I3xx I3yy I3zz
I{1} = [I1xx 0 0; 0 I1yy 0; 0 0 I1zz];
I{2} = [I2xx 0 0; 0 I2yy 0; 0 0 I2zz];
I{3} = [I3xx 0 0; 0 I3yy 0; 0 0 I3zz];

% data
r = cell(1,N);
r0 = cell(1,N);
rcm = cell(1,N);
rcm0 = cell(1,N);
w = cell(1,N);
v = cell(1,N);
vcm = cell(1,N);

syms k

% USER DATA HERE!
r{1} = [q(1)-k; 0; 0];
r{2} = [l(2)*cos(q(2)); l(2)*sin(q(2)); 0];
r{3} = [l(3)*cos(q(2)+q(3)); l(3)*sin(q(2)+q(3)); 0];

%##########################################################################
% moving frame algorithm - planar xy case :::::::::::::::::::::::::::::::::
prev_r  = [0;0;0];
prev_r0 = [0;0;0];
prev_w  = [0;0;0];
prev_v  = [0;0;0];

for i = 1:N
    rcm{i} = subs(r{i}, l(i), dcm(i));
    r0{i} = sum([r{1:i}],2);
    rcm0{i} = prev_r0 + rcm{i};
    w{i} = prev_w + (1-sigma(i))*dq(i)*[0;0;1];
	v{i}    = prev_v + sigma(i)*dq(i)*[1;0;0] + cross(w{i},r{i});
	vcm{i}  = prev_v + sigma(i)*dq(i)*[1;0;0] + cross(w{i},rcm{i});
    
    prev_r0 = r0{i};
    prev_r = r{i};
    prev_w = w{i};
    prev_v = v{i};
end
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%##########################################################################

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

params = [  M(1,1);
            I2zz + I3zz + m(2)*dcm(2)^2 + m(3)*dcm(3)^2 + m(3)*l(2)^2;
            m(3) + dcm(3) + l(2);
            M(3,3);           
            m(2)*dcm(2) + m(3)*l(2);
            m(3)*dcm(3)]
            

[M_a, a] = parametrization(M, params);
[c_a, a] = parametrization(c, params);
[g_a, a] = parametrization(g, params);

Y = simplify(jacobian(M_a*ddq+c_a,a))



