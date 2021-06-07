%	compute POTENTIAL ENERGY
%   @Author: Giuseppe Sensolini ~ March 2020
%
% inputs:   m:          symbolic joint mass vector (Nx1)
%           gravity:    symbolic gravity vector (3x1)
%           rcm0:       symbolic vector from base frame to CoM (cell 1xN of 3x1 vectors)
% outputs:  Ucm:        CoM potential energy vector (Nx1)

function Ucm = potential_energy(m, gravity, rcm0)
    N = length(rcm0);
    Ucm = sym('Ucm',[N,1]);
    for i = 1:N
        Ucm(i) = -m(i) * gravity.' * rcm0{i};
    end
end