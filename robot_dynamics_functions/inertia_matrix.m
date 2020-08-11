%	compute INERTIA MATRIX
%   @Author: Giuseppe Sensolini ~ March 2020
% inputs:   T:          kinetic energy (scalar) 
%           q:          symbolic joint positions vector (Nx1)
%           dq:         symbolic joint velocity vector (Nx1)
%           T:          kinetic energy (scalar)
% outputs:  M:          inertia matrix (Nx1)
function M = inertia_matrix(T, q, dq)
    N = length(q);
    M = sym('M',[N,N]);
    for i = 1:N
        for j = i:N
            M(i,j) = diff( diff(T, dq(i)), dq(j) );
            if i~=j 
                M(j,i) = M(i,j);
            end
        end
    end
    M = simplify(M);
end
