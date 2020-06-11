function ret = fdyn_kinematics(x, u)
    % dynamics for a point mass
    global Slist M_
    ret_qdd = zeros(13, size(x,2));
    ret     = zeros(16, size(x,2));
    T       = @(x)[0 -sin(x(1)) cos(x(1))*cos(x(2));0 cos(x(1)) sin(x(1))*cos(x(2));1 0 sin(x(2))];
    
    for i = 1:size(x,2)
        Mass = [1 0 0;0 1 0;0 0 1]; % for a manipulator, this varies with configuration
        I    = [1 0 0;0 1 0;0 0 1];
        
        % RR    = esul2rotm(x(4:6,i)', 'ZYX');
        
        % friction, in the direction of sliding (tangential) and normal to
        % it
        if (norm(x(14:16,i)) < 0.1)
            vel_dir        = x(14:16,i);
            orthogonal_dir = vel_dir;
        else
            vel_dir        = x(14:16,i)/norm(x(14:16,i));
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
        
        % let level of the surface be z = 0.0
        % distance into the surface level 
        if (x(10,i) <= 0.0)
            dx            = d + x(10,i);
        else
            dx            = d + x(10,i);
        end                
        
        Rot = eul2rotm(x(11:13,i)', 'ZYX');
        tool_length = Rot*[0 0 0.136]'; % end-effector frame of the tool base
        
        f_vec_ee = x(20:22,i);
        
        % singularities in J_
        % q = IKinSpace_modified(Slist, M_, [[Rot;0 0 0],[x(8:10,i);1]], x(1:7,i), 0.001, 0.001);
        % J_ = JacobianSpace(Slist, x(1:7,i));
        
        T_current = [[Rot;0 0 0],[x(8:10,i);1]];
        
        if norm(FKinSpace(M_, Slist, x(1:7,i))-T_current)>0.1
            q = x(1:7,i);
            % q = IKinSpace_modified(Slist, M_, [[Rot;0 0 0],[x(8:10,i);1]], x(1:7,i), 0.001, 0.001);
        else
            q = x(1:7,i);
        end
        
        J_ = JacobianSpace(Slist, q);
        if (rcond(J_*J_') < 0.01)
            J_m      = J_'/(J_*J_' + 0.001*eye(6));
        else
            J_m      = J_'/(J_*J_');
        end
        
        w    = T(x(11:13,i)) * x(17:19,i);
        Rdot = cross(repmat(w,1,3),Rot);
        gdot = [[Rdot;0 0 0],[x(14:16,i);1]];
        g    = [[Rot;0 0 0],[x(8:10,i);1]];
        xd   = gdot/g;
        xd   = [xd(1:3,end);w];
        
        % x(8:10,i)
        % xd
        
        qd   = J_m * xd;
        % xd_  = [x(17:19,i);w];
        % qd   = Jac_kuka() * x();
        
        xdd  = Mass \ (u(1:3,i) + x(20:22,i));
        wdd  = I \ (u(4:6,i) + cross(-tool_length/2, f_vec_ee));   % I is at the base
        
        ret_qdd(:,i)  = [qd;xdd;wdd];
        
        mu   = 0.5;
        nu   = 0.55;
        k    = real((6*E^2*R*x(22,i)/(1-0.9^2)^2)^(1/3));
        
        F_f  = -mu*k*x(16,i) * vel_dir + 0.3*mu*(2*nu-1)/R*(k*x(16,i)*d + x(22,i)*x(16,i));
        
 
        orthogonal_dir = orthogonal_dir(:,1);
        % compute the acceleration vector in the direction othorgonal to
        % the moving direction.
        
        if (norm(orthogonal_dir) > 0.5)
            acc_orth = (dot(xdd(1:3)',orthogonal_dir)/(norm(orthogonal_dir)^2)) * orthogonal_dir;
        else
            acc_orth = 0;
        end
        
        F_normal = 2*1*norm(x(14:16,i)) * norm(acc_orth)/0.2;
        if isnan(F_normal)
           F_normal = 0; 
        end
        
        %  F_normal
        F        = F_f + k*x(16,i)*[0;0;1] + 0*F_normal*orthogonal_dir;
        
        ret(:,i) = [ret_qdd(:,i); F];
        
    end
    
    


end