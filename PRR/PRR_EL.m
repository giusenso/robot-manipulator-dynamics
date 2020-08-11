%%  PRR planar robot - Dynamics
%   @Giuseppe Sensolini, 10 April 2020

disp("### DYNAMICS ########################################################");

%clear all
close all
clear all
clc

%% Robot parameters________________________________________________________

N = 3  % joint number  

q   = sym('q',[N,1]); assume(q,'real');    	% joints symbolic variables
l   = sym('l',[N,1]); assume(l,'real');   	% joints symbolic variables
d   = sym('d',[N,1]); assume(d,'real');   	% distance of the CoM from the joint axes
dq  = sym('dq',[N,1]); assume(dq,'real');   % first derivative of joint varaibles
ddq = sym('ddq',[N,1]); assume(ddq,'real'); % second derivative of joint varaibles
m   = sym('m',[N,1]); assume(m,'real');     % masses of CoM
u   = sym('u',[N,1]); assume(u,'real');     % external forces
r = cell(1,N);  % vector going from frame i-1 to frame i expressed in frame 0
syms g0

% inertia matrices
syms I1xx I1yy I1zz  I2xx I2yy I2zz  I3xx I3yy I3zz
I = cell(1,N);
I{1} = [I1xx 0 0; 0 I1yy 0; 0 0 I1zz];
I{2} = [I2xx 0 0; 0 I2yy 0; 0 0 I2zz];
I{3} = [I3xx 0 0; 0 I3yy 0; 0 0 I3zz];

% ### USER DATA HERE! ###
gravity = [0; 0; -g0];  % gravity vector
sigma = [1 0 0]  % 0 if revolute, 1 if prismatic     
r{1} = [0; 0; q(1)];
r{2} = [0; l(2); 0];
r{3} = [l(3)*cos(q(2))*sin(q(3)); l(3)*cos(q(3)); -l(3)*sin(q(2))*sin(q(3))];
% #######################

disp("=== MOVING FRAME ALGO ===============================================");
[vcm, w, rcm0] = moving_frame_algo(q, dq, r, l, d, sigma);

disp("=== KINETIC ENERGY ==================================================");
[T, Tcm] = kinetic_energy(q, dq, m, I, vcm, w)

disp("=== INERTIA MATRIX ==================================================");
M = inertia_matrix(T, q, dq)
dM = inertia_diff(M, q, dq);

params = [	M(1,1);
           	I{2}(3,3)+m(3)*l(2)^2+m(2)*d(2);
            m(3)*d(3)^2;
            m(3)*d(3)*l(3);
            I{3}(3,3)];
syms a1 a2 a3 a4 a5 real
M22 = a2 + a5 + 2*a4*cos(q(2)) + a3*(cos(q(2))^2+cos(q(2))^2*sin(q(3))^2);
M32 = a5 + a3*(cos(q(3))^2 + cos(q(2))^2*sin(q(3))^2) + a4*cos(q(3));
M33 = a5 + a3*(1 - sin(q(2))^2*sin(q(3))^2)
M_a = [ a1,     0,      0;
        0,      M22,    M32;
        0,      M32,    M33]

disp("=== CORIOLIS AND CENTRIFUGAL ========================================");
[c, christoff, S, X] = centrifugal_coriolis(M_a, q, dq)

disp("=== GRAVITY =========================================================");
[U, Ucm] = potential_energy(q, m, d, gravity, rcm0)
g = simplify(jacobian(U,q).')


%% Parametrization
a = sym('a',[length(params),1])
Y = simplify( jacobian(M_a*ddq + c , a) )











%