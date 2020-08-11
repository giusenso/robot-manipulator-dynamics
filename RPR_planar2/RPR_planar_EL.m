%%  RP planar robot - Dynamics
%   @Giuseppe Sensolini, 07.2020

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
gravity = [g0; 0; 0];   % gravity vector
sigma = [0 1 0]         % 0 if revolute, 1 if prismatic     
r{1} = [0; 0; 0];
r{2} = [q(2)*cos(q(1)); q(2)*sin(q(1)); 0];
r{3} = [l(3)*cos(q(3)); l(3)*sin(q(3)); 0];
% #######################

disp("=== MOVING FRAME ALGO ===============================================");
[vcm, w, rcm0] = moving_frame_algo2(q, dq, r, l, d, sigma);

disp("=== KINETIC ENERGY ==================================================");
[T, Tcm] = kinetic_energy(q, dq, m, I, vcm, w)

disp("=== INERTIA MATRIX ==================================================");
M = inertia_matrix(T, q, dq)
dM = inertia_diff(M, q, dq);

disp("=== CORIOLIS AND CENTRIFUGAL ========================================");
[c, christoff, S, X] = centrifugal_coriolis(M, q, dq)

disp("=== GRAVITY =========================================================");
[U, Ucm] = potential_energy(q, m, d, gravity, rcm0)
g = simplify(jacobian(U,q).')


%% Parametrization
% a = sym('a',[length(params),1])
% Y = simplify( jacobian(M_a*ddq + c , a) )


function [vcm,w,rcm0] = moving_frame_algo2(q, dq, r, l, d, sigma)
    N = length(r);
    
    rcm = cell(1,N);
    rcm0 = cell(1,N);
    w = cell(1,N);
    v = cell(1,N);
    vcm = cell(1,N);
    
    prev_w  = [0;0;0];
    prev_v  = [0;0;0];
    
    rcm{1} = [0; 0; 0];
    rcm{2} = [(q(2)-d(2))*cos(q(1)); (q(2)-d(2))*sin(q(1)); 0];
    rcm{3} = [(l(3)-d(3))*cos(q(3)); (l(3)-d(3))*sin(q(3)); 0];
    
    rcm0{1} = [0; 0; 0];
    rcm0{2} = [(q(2)-d(2))*cos(q(1)); (q(2)-d(2))*sin(q(1)); 0];
    rcm0{3} = r{2}+[(l(3)-d(3))*cos(q(3)); (l(3)-d(3))*sin(q(3)); 0];
    
    for i = 1:N
        w{i}    = prev_w + (1-sigma(i))*dq(i)*[0;0;1];
        v{i}    = prev_v + sigma(i)*dq(i)*(r{i}/q(i)) + cross(w{i},r{i});
        vcm{i}  = prev_v + sigma(i)*dq(i)*(r{i}/q(i)) + cross(w{i},rcm{i});

        prev_w  = w{i};
        prev_v  = v{i};
    end
    print_output(vcm, w, rcm0)
end


function ret = print_output(vcm, w, rcm0)
    for i = 1:length(w)
        fprintf("rcm0%d:\n", i);
        disp(rcm0{i})
        
        fprintf("w%d:\n", i);
        disp(w{i})
        
        fprintf("vcm%d:\n", i);
        disp(vcm{i})
        
        fprintf("________________________________________________\n\n");
    end
end







%