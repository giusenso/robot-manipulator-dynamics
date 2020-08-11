%	compute CENTRIFUGAL and CORIOLIS terms
%   @Author: Giuseppe Sensolini ~ April 2020
% inputs:   M:          inertia matrix (NxN)
%           q:          symbolic joint positions vector (Nx1)
%           dq:         symbolic joint velocity vector (Nx1)
% outputs:  c:          centrifugal and coriolis vector (Nx1)
%           christoff:  christoffel matrices (Nx1 cell containing NxN matrices)
%           S:          skew symmetric factorization (NxN)
%           X:          non-skew symmetric factorization (NxN)
function [c, christoff, S, X] = centrifugal_coriolis(M, q, dq)
    N = length(q);
	c = sym('c',[N,1]); assume(c,'real');
    christoff = cell(1,N);
    S = sym('S',[N,N]);
    for i = 1:1:N
        m_ = jacobian(M(:,i),q);
        christoff{i} = simplify(0.5*(m_ + m_.' - diff(M,q(i))));
        c(i) = dq.'*christoff{i}*dq;   
        for j = 1:N
            S(i,j) = christoff{i}(:,j).' * dq;
        end
    end
    c = simplify(c);
    S = simplify(S);
    X = simplify( jacobian(c,dq) );
end