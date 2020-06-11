function ret = fdyn(x, u)
    % dynamics for a point mass
    ret_qdd = zeros(6, size(x,2));
    ret     = zeros(9, size(x,2));
    
    for i = 1:size(x,2)
        M    = [1 0 0;0 1 0;0 0 1]; % for a manipulator, this varies with configuration
        I    = [1 0 0;0 1 0;0 0 1];
        
        
        %         RR    = eul2rotm(x(4:6,i)', 'ZYX');
        
        % friction, in the direction of sliding (tangential) and normal to
        % it
        if (norm(x(7:9,i)) < 0.1)
            vel_dir        = x(7:9,i);
            orthogonal_dir = vel_dir;
        else
            vel_dir        = x(7:9,i)/norm(x(7:9,i));
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
        u_z           = (dot(u(1:3,i),spring_dir)/(norm(spring_dir)^2)) * spring_dir;
        
        % current manipulator end- effector velocity.
        % curr_vel_spring_dir  = (dot(x(7:9,i),spring_dir)/(norm(spring_dir)^2)) * spring_dir;
        
        
        % surface deformation resulting from u_z, dx in the direction of
        % spring_dir. Quasi static model. 
        d             = (9*norm(u_z)^2/(16*E^2*R))^(1/3);
        
        % let level of the surface be z = 0.0
        % distance into the surface level 
        if (x(3,i) <= 0.0)
            dx            = d + x(3,i);
        else
            dx            = d + x(3,i);
        end                
        
        Rot = eul2rotm(x(4:6,i)', 'ZYX');
        tool_length = Rot*[0 0 0.2]'; % end-effector frame of the tool base
        
        f_vec_ee = x(13:15,i);
        
        xdd  = M \ (u(1:3,i) + x(13:15,i));
        wdd  = I \ (u(4:6,i) + cross(-tool_length, f_vec_ee));   % I is at the base
        
        ret_qdd(:,i)  = [xdd;wdd];
        
        mu = 0.1;
        nu = 0.9;
        k  = real((6*E^2*R*x(15,i)/(1-0.9^2)^2)^(1/3));
        
        F_f      = -mu*k*x(9,i) * vel_dir + 0.3*mu*(2*nu-1)/R*(k*x(9,i)*d + x(15,i)*x(9,i));
        
 
        orthogonal_dir = orthogonal_dir(:,1);
        if (norm(orthogonal_dir) > 0.5)
            acc_orth = (dot(xdd(1:3)',orthogonal_dir)/(norm(orthogonal_dir)^2)) * orthogonal_dir;
        else
            acc_orth = 0;
        end
        
        F_normal = 2*1*norm(x(7:9,i))*norm(acc_orth)/0.2;
        if isnan(F_normal)
           F_normal = 0; 
        end
        
        %  F_normal
        F        = F_f + k*x(9,i)*[0;0;1] + 0*F_normal*orthogonal_dir;
        
        ret(:,i) = [ret_qdd(:,i); F];
        
    end

end