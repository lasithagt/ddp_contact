function dyn = contact_model_dyn(x_r, x_c, u_c)
    n = size(x_c,2);
    dyn = zeros(size(x_c,1),n);
    
    for i = 1:n
        RR    = eul2rotm(x_r(4:6,i)', 'ZYX');

        % friction, in the direction of sliding (tangential) and normal to
        % it
        if (norm(x_r(7:9,i)) == 0)
            vel_dir        = x_r(7:9,i);
            orthogonal_dir = vel_dir;
        else
            vel_dir        = x_r(7:9,i)/norm(x(7:9,i));
            orthogonal_dir = null(vel_dir');
        end

        % friction coefficients 
        mu_fric = 0.1;
        mu_orth = 0.05;

        % material properties of the surface
        E  = 1e1; 
        R  = 10/1000; % in mm

        spring_dir = [0 0 1]';

        % the force component in the direction of spring_dir, projection to
        % the spring_dir
        u_z           = (dot(u_c(1:3,i),spring_dir)/(norm(spring_dir)^2)) * spring_dir;

        % current manipulator end- effector velocity.
        % curr_vel_spring_dir  = (dot(x(7:9,i),spring_dir)/(norm(spring_dir)^2)) * spring_dir;


        % surface deformation resulting from u_z, dx in the direction of
        % spring_dir. Quasi static model. 
        d             = (9*norm(u_z)^2/(16*E^2*R))^(1/3);

        % let level of the surface be z = 0.0
        % distance into the surface level 
        if (x_r(3,i) <= 0.0)
            dx  = d + x_c(3,i);
        else
            dx  = 0;
        end

        % rate at which the surface will move
        dx_dot        = -dx / 10;
        F             = norm(u_z);
        if (F > 0.01)
            F_dot     = (6*E^2*R*(F^2/(E^2*R))^(2/3))/(6^(2/3)*F) * (dx_dot);
        else
            F_dot     = 0;
        end
        dyn(3,i)   = F_dot;


        % force generated from friction.
        % max_normal_fricational_force = ret_fd(3,i)*orthogonal_dir(1:2,1)*mu_orth;

        if (norm(orthogonal_dir(1:3,1))) < 0.1
            normal_fric = [0 0 0]';
        else
            normal_fric   = (dot(RR*spring_dir,orthogonal_dir(1:3,1))/(norm(orthogonal_dir(1:3,1))^2)) * orthogonal_dir(1:3,1);
        end

        dyn(1:2,i) = -F_dot*vel_dir(1:2)*mu_fric - normal_fric(1:2);
    end

end