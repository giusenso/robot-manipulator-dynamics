%%  3R robot
%   Euler-Lagrange Dynamics
%   04.2020, @Giuseppe Sensolini

clc
close all
clear all 

%% Robot parameters _______________________________________________________
disp('::: INPUT DATA :::::::::::::::::::::::::::::::::::');

syms a alpha d theta
sigma = [0 0 0];  % 0 if revolute, 1 if prismatic
N = 3;          % joint number                   
q = sym('q',[N,1]), assume(q,'real');	% joints symbolic variables
l = sym('l',[N,1]), assume(l,'real');	% joints symbolic variables

%% Denavit-Hartemberg table (set up parameters here!) 
disp('::: DH TABLE ::::::::::::::::::::::::::::::::::::::');
% [alpha a d theta]
DHTABLE = [     pi/2    0       l(1)	q(1);
                0       l(2)	0       q(2);
                0       l(3)	0       q(3)]

%% General Denavit-Hartenberg trasformation matrix
DH = [ cos(theta) -sin(theta)*cos(alpha)  sin(theta)*sin(alpha) a*cos(theta);
        sin(theta)  cos(theta)*cos(alpha) -cos(theta)*sin(alpha) a*sin(theta);
          0             sin(alpha)             cos(alpha)            d;
          0               0                      0                   1      ];

%% Direct Kinematics ______________________________________________________
disp("#####################################################################");
disp("### KINEMATICS ######################################################"); 
% Build transformation matrices for each link
A = cell(1,N);
for i = 1:N
    alpha   = DHTABLE(i,1);
    a       = DHTABLE(i,2);
    d       = DHTABLE(i,3);
    theta   = DHTABLE(i,4);
    A{i} = subs(DH);
end
clear alpha a d theta
%% Build base-to-end-effector transformation matrix 
T = eye(4);         % transformation matrix from base to e-e
R = cell(1,N);      % rotation matrices from i-1 to i
r = cell(1,N);      % vector going from link i-1 to i expressed in i-1 frame
rcm0 = cell(1,N);   % vector going from base to the i-th CoM expressed in base frame

% vector going from frame i to the i-th CoM expressed in frame i
syms Am Fm Cm Dm Em real
rcm = [ {[Am;-Fm;0]}, {[-Cm;0;0]}, {[-Dm;0;Em]} ];

for i = 1:N
    T =  simplify(T*A{i});
    R{i} = A{i}(1:3,1:3);
    r{i} = A{i}(1:3,4);
    rcm0{i} = T(1:3,4) + rcm{i};
end

%% DYNAMICS _______________________________________________________________
disp("#####################################################################");
disp("### DYNAMICS ########################################################"); 

syms g0 real;
gravity = [0; 0; -g0];  % gravity vector

dq  = sym('dq',[N,1]); assume(dq,'real');   % first derivative of joint varaibles
ddq = sym('ddq',[N,1]); assume(ddq,'real'); % second derivative of joint varaibles
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

% angular and linear velocity
prev_w = zeros(3,1);
prev_v = zeros(3,1);
w = cell(1,N);
v = cell(1,N);
vcm = cell(1,N);

for i = 1:N
    w_ = simplify( prev_w + (1-sigma(i))*dq(i)*[0;0;1] );  % w{i} wrt frame i-1
    w{i} = simplify( R{i}.' * w_ );
    v{i} = simplify( R{i}.'*(prev_v + sigma(i)*dq(i)*[0;0;1] + cross(w_,r{i})) );
    vcm{i} = simplify( v{i} + sigma(i)*dq(i)*[0;0;1] + cross(w{i},rcm{i}) );
    prev_w = w{i};
    prev_v = v{i};
end

disp("::: KINETIC ENERGY ::::::::::::::::::::::::::::::::::::::::::::::::::");
[T, Tcm] = kinetic_energy(q, dq, m, I, vcm, w)

disp("::: INERTIA MATRIX ::::::::::::::::::::::::::::::::::::::::::::::::::");
M = inertia_matrix(T, q, dq)
dM = inertia_diff(M, q, dq);
%M_a = parametrization(M,[]);

disp("::: CORIOLIS AND CENTRIFUGAL ::::::::::::::::::::::::::::::::::::::::");
[c, christoff, S, X] = centrifugal_coriolis(M, q, dq) 

disp("::: GRAVITY :::::::::::::::::::::::::::::::::::::::::::::::::::::::::");
%[U, Ucm] = potential_energy(q, m, dcm, gravity, rcm0);
%g = simplify(jacobian(U,q).')

disp("::: DYNAMIC MODEL::::::::::::::::::::::::::::::::::::::::::::::::::::");
%simplify( M * ddq + c + g )


