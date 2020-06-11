function qdd = fdyn(parms, x, u)
    
    M = M_kuka(parms, x(1:7));
    C = C_kuka(parms, x(1:7),x(8:end));
    G = G_kuka(parms, x(1:7));
    % M = eye(7);
    size(C)
    
    kv        = diag([4.0, 1.5, 1.0, 0.8, 0.8, 0.3, 0.05]);
    k_static  = (0.0) * diag([0.02 0.02 0.01 0.07 0.01 0.01 0.001]);

    Ff        =   -k_static * sign(x(14:20));
    Fv        =   -kv * x(14:20);


    grav_comp =   -G_kuka(inertial_params, x(1:7))'; 
    qdd       =   kukaAnalytical_qdd(inertial_params,x(1:7), x(14:20), Ff+Fv+grav_comp + u(7:end,i) - J'*u(1:6,i));
        
    qdd = M \ (u - G' - C);
end