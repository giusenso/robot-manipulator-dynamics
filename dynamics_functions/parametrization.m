%   Parametrization
%   substitute groups of constant inside a matrix
%   @Author: Giuseppe Sensolini ~ March 2020
%	inputs:     M:          Matrix or vector
%               params:     Symbolic params to be substituted (Qx1)
%   outputs:    M_a:        Matrix expressed in the new parameters
%               a:          symbolic vector (Qx1)
function [M_a,a] = parametrization(M, params)
    a = sym('a',[length(params),1]);
    M_a = simplify(subs(M, params, a));
end
