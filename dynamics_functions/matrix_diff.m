

function dJ =  matrix_diff(J, q, dq)
    N = size(J,2);
    M = size(J,1);
    for k = 1:M
        for j = 1:N
            for i = 1:M
                if i==1
                    dJ(k,j) = diff(J(k,j),q(i)) * dq(i);
                else
                    dJ(k,j) = dJ(k,j) + diff(J(k,j),q(i)) * dq(i);
                end
                dJ(k,j) = simplify(dJ(k,j));
            end
        end
    end
end