%   compute KINETIC ENERGY
%   @Author: Giuseppe Sensolini ~ March 2020
%
% inputs:   q:          symbolic joint positions vector (Nx1)
%           dq:         symbolic joint velocity vector (Nx1)
%           m:          symbolic joint mass vector (Nx1)
%           I:          symbolic joint inertia vector (Nx1)
%           vcm:        symbolic CoM velocity (1xN cell containing 3x1 vectors)
%           w:          symbolic angular velocity (1xN cell containing 3x1 vectors)
% outputs:  T:          total manipulator kinetic energy (scalar)
%           Tcm:        CoM kinetic energy vector (Nx1)

function [T, Tcm] = kinetic_energy(q, dq, m, I, vcm, w)
    N = length(q);
    Tcm = sym('Tcm',[N,1]); assume(Tcm,'real');
    for i = 1:N
        transl_energy = 0.5 * m(i) * (vcm{i}.'*vcm{i});
        rot_energy = 0.5 * (w{i}.' * I{i} * w{i});
        Tcm(i) = simplify( transl_energy + rot_energy );
    end
    T = simplify( sum(Tcm) );
end