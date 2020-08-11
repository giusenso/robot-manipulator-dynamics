%%  RP planar robot - Dynamics
%   @Giuseppe Sensolini, 10 April 2020

disp("### DYNAMICS ########################################################");

%clear all
close all
clear all
clc

%% Robot parameters________________________________________________________

N = 2  % joint number  

q   = sym('q',[N,1]); assume(q,'real');    	% joints symbolic variables
l   = sym('l',[N,1]); assume(l,'real');   	% joints symbolic variables
d   = sym('d',[N,1]); assume(d,'real');   	% distance of the CoM from the joint axes
dq  = sym('dq',[N,1]); assume(dq,'real');   % first derivative of joint varaibles
ddq = sym('ddq',[N,1]); assume(ddq,'real'); % second derivative of joint varaibles
m   = sym('m',[N,1]); assume(m,'real');     % masses of CoM
u   = sym('u',[N,1]); assume(u,'real');     % external forces
r = cell(1,N);	% vector going from frame i-1 to frame i expressed in frame 0
syms g0

% inertia matrices
syms I1xx I1yy I1zz  I2xx I2yy I2zz
I = cell(1,N);
I{1} = [I1xx 0 0; 0 I1yy 0; 0 0 I1zz];
I{2} = [I2xx 0 0; 0 I2yy 0; 0 0 I2zz];

% ### USER DATA HERE! ###
gravity = [0; -g0; 0];  % gravity vector
sigma = [0 1]           % 0 if revolute, 1 if prismatic     
r{1} = [l(1)*cos(q(1)); l(1)*sin(q(1)); 0];
r{2} = [-q(2)*sin(q(1)); q(2)*cos(q(1)); 0];
% #######################
%%

disp("=== MOVING FRAME ALGO ===============================================");
[vcm, w, rcm0] = moving_frame_algo_v11(q, dq, r, l, d, sigma);

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



%%
% for this specific problem we need to sligthly modify the moving frame
% algorithm becouse of the particular choice of d1 and d2
% (look at the robot scheme in the attached figure)
function [vcm,w,rcm0] = moving_frame_algo_v11(q, dq, r, l, d, sigma)
    N = length(r);
    
    r0 = cell(1,N);
    rcm = cell(1,N);
    rcm0 = cell(1,N);
    w = cell(1,N);
    v = cell(1,N);
    vcm = cell(1,N);
    
    prev_r0 = [0;0;0];
    prev_w  = [0;0;0];
    prev_v  = [0;0;0];
    
    for i = 1:N
        r0{i}   = sum([r{1:i}],2);
        if sigma(i)==0
            rcm{i}  = subs(r{i}, l(i), d(i));
        else
            rcm{i}  = subs(r{i}, q(i), q(i)-d(i));
        end
        rcm0{i} = prev_r0 + rcm{i};
        w{i}    = prev_w + (1-sigma(i))*dq(i)*[0;0;1];
        v{i}    = prev_v + sigma(i)*dq(i)*(r{i}/q(i)) + cross(w{i},r{i});
        vcm{i}  = prev_v + sigma(i)*dq(i)*(r{i}/q(i)) + cross(w{i},rcm{i});

        prev_r0 = r0{i};
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