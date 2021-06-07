% PRINT UTILITY

function print_dynamics(sigma, vcm, w, rcm0, Tcm, M, c, christoff, S, Ucm, g)
    fprintf('  ~ Robot Manipulator Dynamics ~\n\n  Developed by Giuseppe Sensolini (github: @giusenso)\n  May 2020\n');
    fprintf('  Brief: algorithms for computing Euler-Lagrange dynamics\n\n\n')
    
    fprintf('robot: %s\n\n', joint_types(sigma));

    disp("=== MOVING FRAME ================================");
    N = length(sigma);
    for i = 1:N
        fprintf("rcm0%d =\n", i); disp(rcm0{i})     
        fprintf("w%d =\n", i); disp(w{i})
        fprintf("vcm%d =\n", i); disp(vcm{i})
        fprintf("Tcm(%d) = ", i); disp(Tcm(1));
        if i<N
            fprintf("__________________________________\n\n");
        end
    end
    
    disp("=== KINETIC ENERGY ==============================");
    fprintf("\nT = "); disp(sum(Tcm));

    disp("=== INERTIA MATRIX ==============================");
    fprintf("\nM(q) = \n"); disp(M);
    
    disp("=== Christhoffel symbols ========================");
    for i=1:N
        disp(['C(',num2str(i),') = ']); disp(christoff{i});
    end
    
    disp("=== CORIOLIS and CENTRIFUGAL ====================");
    fprintf("\nc(q,dq) = \n"); disp(c);
    fprintf("S(q,dq) = \n"); disp(S);
    
    disp("=== POTENTIAL ENERGY ============================");
    fprintf("\nU = "); disp(sum(Ucm));
    
    disp("=== GRAVITY =====================================");
    fprintf("\ng(q) = \n"); disp(g);
    
    disp("=== ROBOT DYNAMICS ==============================");
    fprintf('\n');
    u = sym('u',[N,1]); assume(u,'real');       % external forces
    ddq = sym('ddq',[N,1]); assume(ddq,'real'); % acc inputs
    dyn = simplify(M*ddq + c + g);
    disp(dyn==u) 
end


function robot = joint_types(sigma)
    robot = '';
    for i=1:length(sigma)
       if sigma(i)==0
           robot(i)='R';
       else
           robot(i)='P';
       end
    end
end