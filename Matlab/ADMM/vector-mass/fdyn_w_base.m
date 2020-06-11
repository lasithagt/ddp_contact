function ret = fdyn_w_base(x, u)
    % dynamics for a point mass
    ret_qdd = zeros(6, size(x,2));
    ret     = zeros(12, size(x,2));
            
    T = @(x)[0 -sin(x(1)) cos(x(1))*cos(x(2));0 cos(x(1)) sin(x(1))*cos(x(2));1 0  sin(x(2))];
    T_dot = @(x,y)[0 -cos(x(1))*y(1) -cos(x(1))*sin(x(2))*y(2)-sin(x(1))*cos(x(2))*y(1);...
        0 -sin(x(1))*y(1) -sin(x(1))*sin(x(2))*y(2)+cos(x(1))*cos(x(2))*y(1);...
        0 0 cos(x(2))*y(2)];
    
    for i = 1:size(x,2)
        M    = [1 0 0;0 1 0;0 0 1]; % for a manipulator, this varies with configuration
        I    = [01 0 0;0 01 0;0 0 1];
        
        spring_dir = [0 0 1]';
        u_z        = (dot(u(1:3,i),spring_dir)/(norm(spring_dir)^2)) * spring_dir;

        R = eul2rotm(x(4:6,i)', 'ZYX');
        
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
        
        % current manipulator end- effector velocity.
        % curr_vel_spring_dir  = (dot(x(7:9,i),spring_dir)/(norm(spring_dir)^2)) * spring_dir;
        
        r_tool   = [0 0 0.136]'/2; % end-effector frame of the tool base
        
        f_vec_ee = x(19:21,i);
        xdd      = M \ (u(1:3,i) - f_vec_ee); % this is right
        I_       = R*I*R';
        w        = T(x(4:6,i)) * x(13:15,i);
        wdd      = I_ \ (u(4:6,i) + cross(R*r_tool, -f_vec_ee) - 0*cross(w,I_*w));   % this is right
        
        
        % project the dynamics of the CoM to the base (the base movement).
        % To penalize the movement of the base
        % w        = T(x(4:6,i)) * x(13:15,i);
        % w_dot    = T_dot(x(4:6,i),x(13:15,i)) * x(13:15,i) + T(x(4:6,i)) * wdd;
        eul_dd = T(x(4:6,i)) \ (wdd - T_dot(x(4:6,i),x(13:15,i)) * x(13:15,i));
        x_b_dot  = x(10:12,i) + cross(w,-R*r_tool);
        % x_b_dot  = x(16:18,i);% + cross(w,-R*r_tool);
        x_b_ddot = xdd(1:3) + cross(wdd,-R*r_tool) + cross(w, cross(w,-R*r_tool)) + 2*cross(w,x_b_dot);

        % project the dynamics of the CoM to the ee .
        % To penalize the movement of the base
%         w        = T(x(4:6,i)) * x(13:15,i);
%         w_dot_e  = T_dot(x(4:6,i),x(13:15,i)) * x(13:15,i) + T(x(4:6,i)) * wdd;
%         % x_e_dot  = x(10:12,i) + cross(w,R*r_tool);
%         x_e_dot  = x(16:18,i) + cross(w,R*r_tool);
%         x_e_ddot = xdd(1:3) + cross(w_dot_e,R*r_tool) + cross(w, cross(w,R*r_tool)) + 2*cross(w,x_e_dot);
%            
        % ret_qdd(:,i)  = [x_e_ddot;w_dot_e];
        ret_qdd(:,i)  = [xdd;eul_dd];
        
        % material properties of the surface
        E       = 500; 
        Rd      = 10/100; % in mm
        s_n     = [0 0 1]'; % for time being 
        mu      = 0.1;
        nu      = 0.9;
        
        d    = (9*norm(f_vec_ee)^2/(16*E^2*Rd))^(1/3);
        k    = real((6*E^2*Rd*x(21,i)/(1-0.9^2)^2)^(1/3)); 
        
        F_f  = -mu*k*x(12,i) * vel_dir + 0.3*mu*(2*nu-1)/Rd*(k*x(12,i)*d + x(21,i)*x(12,i));
        

        
        F        = F_f + k*x(12,i)*[0;0;1];
        
        ret(:,i) = [ret_qdd(:,i); x_b_ddot; F];
        
    end

end