function ret = fdyn_dynamics_ext(x, u)
    % dynamics for a point mass
    global inertial_params
    ret_qdd = zeros(13, size(x,2));
    ret     = zeros(16, size(x,2));
    T       = @(x)[0 -sin(x(1)) cos(x(1))*cos(x(2));0 cos(x(1)) sin(x(1))*cos(x(2));1 0 sin(x(2))];
    
    for i = 1:size(x,2)
        Mass = [1 0 0;0 1 0;0 0 1]; % for a manipulator, this varies with configuration
        I    = [1 0 0;0 1 0;0 0 1];
        
        
        % friction, in the direction of sliding (tangential) and normal to
        if (norm(x(21:23,i)) < 0.1)
            vel_dir        = x(21:23,i);
            orthogonal_dir = vel_dir;
        else
            vel_dir        = x(21:23,i)/norm(x(21:23,i));
            if isnan(vel_dir(1)) || isinf(vel_dir(1))
                orthogonal_dir = [0 0 0]';
            else
                orthogonal_dir = null(vel_dir');
            end
        end
        
        
        % material properties of the surface
        E  = 100; 
        R  = 10/100; % in mm
                
        spring_dir = [0 0 1]';
        
        % the force component in the direction of spring_dir, projection to
        % the spring_dir
        u_z        = (dot(u(1:3,i),spring_dir)/(norm(spring_dir)^2)) * spring_dir;
        
        % current manipulator end- effector velocity.
        % curr_vel_spring_dir  = (dot(x(7:9,i),spring_dir)/(norm(spring_dir)^2)) * spring_dir;
        
        
        % surface deformation resulting from u_z, dx in the direction of
        % spring_dir. Quasi static model. 
        d             = (9*norm(u_z)^2/(16*E^2*R))^(1/3);
        
        Rot = eul2rotm(x(11:13,i)', 'ZYX');
        tool_length = Rot*[0 0 0.136]'; % end-effector frame of the tool base
        
        f_vec_ee = x(27:29,i);
                
        q = x(1:7,i);

        R_        = eul2rotm([x(11,i),x(12,i),x(13,i)], 'ZYX');
        T        = [[R_ ;[0 0 0]], [x(8,i);x(9,i);x(10,i);1]];
        
        J         = Jac_kuka(q);
        
        kv            = diag([4.0, 0.5, 0.5, 0.5, 0.5, 0.3, 0.1]);
        k_static      = (0.0) * diag([0.02 0.02 0.01 0.07 0.01 0.01 0.001]);
        
        Ff            =   -k_static * sign(x(14:20,i));
        Fv            =   -kv * x(14:20,i);
        % grav_comp     =   gravityTorque(lbr, x(1:7,i))
        % qdd           =   forwardDynamics(lbr, x(1:7,i), x(14:20,i), Ff+grav_comp+Fv+u(7:end,i));
        
        grav_comp     =   -G_kuka(inertial_params, x(1:7,i))'; 
        qdd           =   kukaAnalytical_qdd(inertial_params,x(1:7,i), x(14:20,i), Ff+Fv+grav_comp + u(7:end,i) - J'*u(1:6,i));
        %         M_k   = eye(7); % just for testing
        %         qdd   = M_k \ (u(7:end,i)+J_'*u(1:6,i)*0);
        
        xdd  = Mass \ (u(1:3,i) + x(27:29,i));
        wdd  = I \ (u(4:6,i) + cross(-tool_length/2, f_vec_ee));   % I is at the base
        
        ret_qdd(:,i)  = [qdd;xdd;wdd];
        
        mu   = 0.5;
        nu   = 0.55;
        k    = real((6*E^2*R*x(29,i)/(1-0.9^2)^2)^(1/3));
        
        F_f  = -mu*k*x(23,i) * vel_dir + 0.3*mu*(2*nu-1)/R*(k*x(23,i)*d + x(29,i)*x(23,i));
        
 
        orthogonal_dir = orthogonal_dir(:,1);
        % compute the acceleration vector in the direction othorgonal to
        % the moving direction.
        
        if (norm(orthogonal_dir) > 0.5)
            acc_orth = (dot(xdd(1:3)',orthogonal_dir)/(norm(orthogonal_dir)^2)) * orthogonal_dir;
        else
            acc_orth = 0;
        end
        
        F_normal = 2*1*norm(x(21:23,i)) * norm(acc_orth)/0.2;
        if isnan(F_normal)
           F_normal = 0; 
        end
        
        %  F_normal
        F        = F_f + k*x(23,i)*[0;0;1] + 0*F_normal*orthogonal_dir;
        ret(:,i) = [ret_qdd(:,i); F];
        
    end
    
    


end