function [q_vec] = fk_test(x,q_0)

    T_el       = @(x) [0 -sin(x(1)) cos(x(1))*sin(x(2));0 cos(x(1)) sin(x(1))*sin(x(2));1 0 cos(x(2))];

    R            = eul2rotm([x(4,1),x(5,1),x(6,1)], 'ZYZ');
    T            = [[R ;[0 0 0]], [x(1,1);x(2,1);x(3,1);1]];
    eomg = 0.01; 
    ev   = 0.001;
    [Slist, M_K] = manipulator_POE();
    [q, s]       = IKinSpace_modified(Slist, M_K, T, q_0, eomg, ev);
    q_0          = q;
    q_vec(:,1)   = q_0;


    for i = 1:size(x,2)

        J        = JacobianSpace(Slist, q);

        J_       = J;
        
        % singularities in J_
        if (rcond(J_*J_') < 0.00001)
            J_m      = J_'/(J_*J_' + 0.001*eye(6));
        else
            J_m      = J_'/(J_*J_');
        end
        
        x_d      = x(7:12,i);
        x_d(4:6) = T_el(x(4:6,i)) * x_d(4:6);
        qd       = J_m  * x_d;

        % q        = q + qd * 0.01;
        
        Tsb_     = FKinSpace(M_K, Slist, q);
        
        R        = eul2rotm([x(4,i),x(5,i),x(6,i)], 'ZYZ');
        T        = [[R ;[0 0 0]], [x(1,i);x(2,i);x(3,i);1]];
    
        Vs       = Adjoint(Tsb_) * se3ToVec(MatrixLog6(TransInv(Tsb_) * T));
        
        
        
        q        = q + 1*J_m*Vs;
        
        q_vec(:,i+1) = q;

    end
end