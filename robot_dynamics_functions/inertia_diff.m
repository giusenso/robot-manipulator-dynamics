

function dM = inertia_diff(M, q, dq)
    N = length(q);
    for k = 1:N
        for j = 1:N
            for i = 1:N
                if i==1
                    dM(k,j) = diff(M(k,j),q(i)) * dq(i);
                else
                    dM(k,j) = dM(k,j) + diff(M(k,j),q(i)) * dq(i);
                end
                dM(k,j) = simplify(dM(k,j));
            end
        end
    end
end