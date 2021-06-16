# Robot Manipulator Dynamics
This repository contains Matlab scripts used to compute Euler-Lagrange dynamic models of a large number of robot manipulators.

Euler-Lagrange Dynamics:  ```M(q)*ddq + c(q,dq) + g(q) = u```

**_dynamics_functions_** contains all the functions needed to compute and display the Euler Lagrange dynamics.
In particular:
```
- moving_frame.m
- kinetic_energy.m
- inertia_matrix.m
- centrifugal_coriolis.m
- factorization.m
- potential_energy.m
- print_dynamics.m
```

Make sure to setup the Matlab path correctly in order to use this library, or alternatively copy and paste it at the end of your scripts.
