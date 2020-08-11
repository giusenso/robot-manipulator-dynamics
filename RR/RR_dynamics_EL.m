%%  2R robot
%   Euler-Lagrange Dynamics
%   04.2020, @Giuseppe Sensolini

clc
close all
clear all 

%% Robot parameters _______________________________________________________
disp('::: INPUT DATA :::::::::::::::::::::::::::::::::::');

syms a alpha d theta
sigma = [0 0];  % 0 if revolute, 1 if prismatic
N = 2;              % joint number                   
q = sym('q',[N,1]), assume(q,'real');	% joints symbolic variables
l = sym('l',[N,1]), assume(l,'real');	% joints symbolic variables

%% Denavit-Hartemberg table (set up parameters here!) 
disp("#####################################################################");
disp("### KINEMATICS ######################################################"); 
disp('::: DH TABLE ::::::::::::::::::::::::::::::::::::::');
% [alpha a d theta]
DHTABLE = [     pi/2    0       l(1)    q(1);
                0       l(2)    0       q(2)]

%% General Denavit-Hartenberg trasformation matrix
DH = [ cos(theta) -sin(theta)*cos(alpha)  sin(theta)*sin(alpha) a*cos(theta);
        sin(theta)  cos(theta)*cos(alpha) -cos(theta)*sin(alpha) a*sin(theta);
          0             sin(alpha)             cos(alpha)            d;
          0               0                      0                   1      ];

%% Direct Kinematics ______________________________________________________
% Build transformation matrices for each link
for i = 1:N
    alpha   = DHTABLE(i,1);
    a       = DHTABLE(i,2);
    d       = DHTABLE(i,3);
    theta   = DHTABLE(i,4);
    A{i} = subs(DH);
end
clear alpha a d theta

% Build base-to-end-effector transformation matrix 
T = eye(4);
dcm = sym('dcm',[N,1]); assume(dcm,'real'); % discance of the CoM from the joint axes
R = cell(1,N);
r = cell(1,N);
rcm = cell(1,N);
rcm0 = cell(1,N);
for i = 1:N
    T =  T * A{i};
    T = simplify(T);
    R{i} = A{i}(1:3,1:3);  % rotation matrix from link i-1 to i
    r{i} = A{i}(1:3,4);    % position vector from link i-1 to i
    rcm{i} = subs(r{i}, l(i), dcm(i)); % position vector from link i-1 to i-th CoM 
    rcm0{i} = subs(T(1:3,4),l(i),dcm(i));
end
T

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
syms I1xx I1yy I1zz  I2xx I2yy I2zz  I3xx I3yy I3zz  I4xx I4yy I4zz
I{1} = [I1xx 0 0; 0 I1yy 0; 0 0 I1zz];
I{2} = [I2xx 0 0; 0 I2yy 0; 0 0 I2zz];
I{3} = [I3xx 0 0; 0 I3yy 0; 0 0 I3zz];
I{4} = [I4xx 0 0; 0 I4yy 0; 0 0 I4zz];

% angular and linear velocity
prev_w = zeros(3,1);
prev_v = zeros(3,1);
w = cell(1,N);
v = cell(1,N);
vcm = cell(1,N);

% vector going from frame i to the i-th CoM expressed in frame i
rcm = [ {[0; -l(1)+dcm(1); 0]}, {[-l(2)+dcm(2); 0; 0]} ];

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

params = [  m(2)*dcm(2)^2+I2yy;
            I2xx;
            I1yy;
            I2zz+m(2)*dcm(2)^2]
a = sym('a',[length(params),1]);
M_a = [a(1)*cos(q(2))^2 + a(2)*sin(q(2)^2) + a(3), 0; 0, a(4)]

        
disp("::: CORIOLIS AND CENTRIFUGAL ::::::::::::::::::::::::::::::::::::::::");
[c, christoff, S, X] = centrifugal_coriolis(M_a, q, dq) 

disp("::: GRAVITY :::::::::::::::::::::::::::::::::::::::::::::::::::::::::");
[U, Ucm] = potential_energy(q, m, dcm, gravity, rcm0);
g = simplify(jacobian(U,q).')





%% PLOT
link = [1; 1];
q0 = [0; 0];
dh = subs(DHTABLE, l, link);
dh = double(subs(dh, q, q0));

% [THETA D A ALPHA TYPE]
for i = 1:N
    L(i) = Link([dh(i,4),dh(i,3),dh(i,2),dh(i,1), sigma(i)], 'standard');
end

rob = SerialLink(L)
rob.name = '2R robot';
rob.plot(q0.');
