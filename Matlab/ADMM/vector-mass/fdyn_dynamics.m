function ret = fdyn_dynamics(x, u, RC, K_D)
    
    % dynamics for a point mass
    global inertial_params tool
    ret_qdd = zeros(13, size(x,2));
    ret     = zeros(16, size(x,2));
    
    T_eu  = @(x)[0 -sin(x(1)) cos(x(1))*cos(x(2));0 cos(x(1)) sin(x(1))*cos(x(2));1 0 sin(x(2))];
    T_dot = @(x,y)[0 -cos(x(1))*y(1) -cos(x(1))*sin(x(2))*y(2)-sin(x(1))*cos(x(2))*y(1);...
        0 -sin(x(1))*y(1) -sin(x(1))*sin(x(2))*y(2)+cos(x(1))*cos(x(2))*y(1);...
        0 0 cos(x(2))*y(2)];
    
    for i = 1:size(x,2)
        Mass = [0.3 0 0;0 0.3 0;0 0 0.3]; % for a manipulator, this varies with configuration
        I    = [1 0 0;0 1 0;0 0 1];
        
        
        % friction, in the direction of sliding (tangential) and normal to
        if (norm(x(21:23,i)) < 0.1)
            vel_dir        = x(21:23,i);
        else
            vel_dir        = x(21:23,i)/norm(x(21:23,i));
        end
        
        
        % material properties of the surface
        E  = 100000; 
        R  = 0.0013/2; % in mm
        m  = 0.3;
                
        
        
        % the force component in the direction of spring_dir, projection to
        % the spring_dir
        spring_dir = [0 0 1]';
        u_z        = (dot(u(1:3,i),spring_dir)/(norm(spring_dir)^2)) * spring_dir;
        
        
        % surface deformation resulting from u_z, dx in the direction of
        % spring_dir. Quasi static model. 
        d        = (9*x(29,i)^2/(16*E^2*R))^(1/3);
        
        % Rot = eul2rotm(x(11:13,i)', 'ZYX');
        % tool_length = Rot*[0 0 0.136]'; % end-effector frame of the tool base
        
        f_vec_ee = x(27:29,i);
                
        q         = x(1:7,i);
        R_        = eul2rotm([x(11,i),x(12,i),x(13,i)], 'ZYX');
        T         = [[R_ ;[0 0 0]], [x(8,i);x(9,i);x(10,i);1]];
        
        J         = Jac_kuka(q); % jacobian at the base of the manipulator
        
        kv        = diag([4.0, 1.5, 1.0, 0.8, 0.8, 0.3, 0.05]);
        k_static  = (0.0) * diag([0.02 0.02 0.01 0.07 0.01 0.01 0.001]);
        
        Ff        =   -k_static * sign(x(14:20,i));
        Fv        =   -kv * x(14:20,i);
   
        
        grav_comp =   -G_kuka(inertial_params, x(1:7,i))'; 
        qdd       =   kukaAnalytical_qdd(inertial_params,x(1:7,i), x(14:20,i), Ff+Fv+grav_comp + u(7:end,i) - J'*u(1:6,i));
        
        %         xdd  = Mass \ (u(1:3,i) - x(27:29,i));
        %         wdd  = I \ (u(4:6,i) + cross(tool_length/2, -f_vec_ee));   % I is at the base

        xdd      = Mass \ (-u(1:3,i) + f_vec_ee) ; % this is right
        I_       = R_*I*R_';
        w        = T_eu(x(11:13,i)) * x(24:26,i);
        wdd      = I_ \ (u(4:6,i) + cross(R_*tool/2, -f_vec_ee));   % this is right

        eul_dd   = T_eu(x(11:13,i)) \ (wdd - T_dot(x(11:13,i),x(24:26,i)) * x(24:26,i)); % convert to euler acceleration
        
        % project the centroid dynamics to the end-point
        xdd_e    = xdd + cross(wdd, R_*tool/2) + cross(w, cross(w, R_*tool/2)); % dynamics of the end point
        
        ret_qdd(:,i)  = [qdd;xdd_e;eul_dd];
        
        mu       = 0.4340;
        nu       = 0.55;
        

        k        = (6*E^2*R*abs(x(29,i)))^(1/3);
        
        if (k < 10)
            k = (6*E^2*R*abs(0.01))^(1/3);
        end
        
        % frictional force
        F_f      = mu*k*x(23,i)*vel_dir + 3*mu*(2*nu-1)*(k*x(23,i)*d+x(29,i)*x(23,i))*vel_dir/(10*R) + 14.02*xdd_e;
        
        % orthogonal_dir = orthogonal_dir(:,1);
        %         % compute the acceleration vector in the direction othorgonal to
        %         % the moving direction.
         
        if (norm(K_D(:,i)) == 0)
            K_DIR = K_D(:,i);
        else
            K_DIR    = K_D(:,i) ./ norm(K_D(:,i));
        end
        
        % acc_orth = (dot(xdd(1:3)',K_DIR')) * K_DIR;
        % acc_orth = (dot(xdd(1:3)',orthogonal_dir)) * orthogonal_dir;
        % F_normal_dot = 2*0.3*norm(x(21:23,i)) * norm(acc_orth)/RC(i);
        
        F_normal_dot = 2 * m * x(21:23,i) .* xdd(1:3) / RC(i);
        %         if isnan(F_normal_dot)
        %            F_normal_dot = 0; 
        %         end
     
%         F        = F_f - k*x(23,i)*[0;0;1] - 80*xdd_e(3)*[0;0;1] + 0*F_normal_dot.*K_DIR;
        F        = 0*F_f - k*x(23,i)*[0;0;1] - 0*xdd_e(3)*[0;0;1] + 0*F_normal_dot.*K_DIR;
        ret(:,i) = [ret_qdd(:,i); F];
        
    end
    
    


end