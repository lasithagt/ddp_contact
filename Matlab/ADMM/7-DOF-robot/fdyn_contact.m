function qdd = fdyn_contact(parms, x, u)
    
    qdd = zeros(10,size(x,2));
    ret_fd = zeros(3,size(x,2));
    
    for i=1:size(x,2)
        M = M_kuka(parms, x(1:7,:));
        C = C_kuka(parms, x(1:7),x(8:end));
            
        fk    = FK_kuka(x(1:7,i));
        R     = fk(1:3,1:3);
        
        % friction, in the direction of sliding (tangential) and normal to
        % it
        if (norm(x(7:9,i)) == 0)
            vel_dir        = x(7:9,i);
            orthogonal_dir = vel_dir;
        else
            vel_dir        = x(7:9,i)/norm(x(7:9,i));
            orthogonal_dir = null(vel_dir');
        end
        
        % friction coefficients 
        mu_fric = 0.1;
        mu_orth = 0.05;
        
        % material properties of the surface
        E   = 1e1; 
        Rd  = 10/1000; % in mm
                
        spring_dir = [0 0 1]';
        
        % the force component in the direction of spring_dir, projection to
        % the spring_dir
        
        u_z = (dot(x(15:17,i),spring_dir)/(norm(spring_dir)^2)) * spring_dir;
        
        % current manipulator end- effector velocity.
        % curr_vel_spring_dir  = (dot(x(7:9,i),spring_dir)/(norm(spring_dir)^2)) * spring_dir;
        
        
        % surface deformation resulting from u_z, dx in the direction of
        % spring_dir. Quasi static model. 
        d  = (9*norm(1)^2/(16*E^2*Rd))^(1/3);
        
        % let level of the surface be z = 0.0
        % distance into the surface level 
        if (fk(3,end) >= 0.9)
            dx            = d - (fk(3,end)-0.9);
        else
            dx            = 0;
        end
        
        % rate at which the surface will move
        dx_dot        = dx / 100;
        F             = norm(x(15:17,i));
        if (F >= 0.0)
            F_dot     = (6*E^2*Rd*(F^2/(E^2*Rd))^(2/3))/(6^(2/3)*F) * (dx_dot);
            if isnan(F_dot)
                F_dot = 1;
            end
        else
            F_dot     = 0;
        end
        ret_fd(3,i)   = F_dot;
%         F_dot
%         x(17,i)
        
        % force generated from friction.
        % max_normal_fricational_force = ret_fd(3,i)*orthogonal_dir(1:2,1)*mu_orth;
        
        %         if (norm(orthogonal_dir(1:3,1))) < 0.1
        %             normal_fric = [0 0 0]';
        %         else
        %             normal_fric = (dot(R*spring_dir,orthogonal_dir(1:3,1))/(norm(orthogonal_dir(1:3,1))^2)) * orthogonal_dir(1:3,1);
        %         end
        
        % ret_fd(1:2,i) = -ret_fd(3,i)*vel_dir(1:2)*mu_fric - normal_fric(1:2);
        
        J = Jac_kuka(x(1:7,i));
        J = J(1:3,:);
        
        qdd(:,i) = [M \ (u(:,i) - 0*J'*x(15:17,i));0;0;ret_fd(3,i)];
    end
end