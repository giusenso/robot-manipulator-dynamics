%	compute POTENTIAL ENERGY
%   @Author: Giuseppe Sensolini ~ March 2020
% inputs:   q:          symbolic joint positions vector (Nx1)
%           m:          symbolic joint mass vector (Nx1)
%           dcm:        symbolic CoM distance (from previous joint) vector (Nx1)
%           gravity:    symbolic gravity vector (3x1)
%           rcm0:       symbolic vector from base frame to CoM (cell 1xN of 3x1 vectors)
% outputs:  U:          total manipulator potential energy (scalar)
%           Ucm:        CoM potential energy vector (Nx1)
function [U, Ucm] = potential_energy(q, m, dcm, gravity, rcm0)
    N = length(q);
    U = 0;
    Ucm = sym('Ucm',[N,1]);
    for i = 1:N
        Ucm(i) = -m(i) * gravity.' * rcm0{i};
        U = U + Ucm(i);
    end
    U = simplify(U);
end