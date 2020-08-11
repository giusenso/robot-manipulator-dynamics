%   MOVING FRAME ALGORITHM V3 for Planar robots
%   Compute linear velocity, angular velocity and position of all CoM of
%   every planar robot placed on xy plane.
%   Assumption: CoM are placed along the link axis.
%   @Author: Giuseppe Sensolini ~ April 2020
% 
% inputs:   q:          symbolic joint positions vector (Nx1)
%           dq:         symbolic joint velocity vector (Nx1)
%           r:          cell (1xN) containing vectors (3x1) going from frame i-i to i
%                       expressed in frame i-1
%           l:          symbolic link length vector (Nx1)
%           d:          symbolic distance from frame i-1 to the i-th CoM
%           sigma:      boolean joint type vector. 0=rev, 1=prism (1xN)
% outputs:  vcm:        cell (1xN) containing velocity vectors (3x1) of every CoM
%           w:          cell (1xN) containing angular velocity vectors (3x1) of every CoM
%           rcm0:       cell (1xN) containing vectors (3x1) going from frame 0 to
%                       the i-th CoM

function [vcm,w,rcm0] = moving_frame_algo_v3(q, dq, r, l, d, sigma)
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
            rcm{i}  = subs(r{i}, l(i), l(i)-d(i));
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




